function obj = Calculate_HeatTransfer_Error(obj)
% CALCULATE_HEATTRANSFER_ERROR computes the error in the heat transfer

%   Computation with the inputs in the obj.HeatTransfer_params.Error
%   struct. The obj.CalculateHeatTransferErrorMethod attribute controls
%   the computation method. There are two methods available:
%   - Moffat: based on the process described in [1].
%   - Montecarlo: based on the process described in [2]

%   See also the routines in @PIRT

%   References:
%   [1] Moffat R J (1988) Describing the uncertainties in experimental results 
%       Exp. Therm. Fluid Sci. 1 3–17 10.1016/0894-1777(88)90043-X
%   [2] Minkina, W., & Dudzik, S. (2009). Infrared thermography: 
%       errors and uncertainties. John Wiley & Sons.

%   Author: I. Robledo, J. Alfaro, R. Castellanos
%   Copyright 2023 Universidad Carlos III de Madrid.


if strcmp(obj.CalculateHeatTransferErrorMethod,'Moffat')
    obj = moffat_method(obj);
end

if strcmp(obj.CalculateHeatTransferErrorMethod,'Montecarlo')
    obj = montecarlo_method(obj);
end

end

function obj = moffat_method(obj)
%% Obtain data

if ~isempty(obj.result)
    if isfield(obj.result,'Thot_new')
        Thot = obj.result.Thot_new;
    else
        Thot = obj.Thot;
    end
    if isfield(obj.result,'Tcold_new')
        Tcold = obj.result.Tcold_new;
    else
        Tcold = obj.Tcold;
    end
else
    Thot = obj.Thot;
    Tcold = obj.Tcold;
end

s               = obj.HeatTransfer_params.HFS.s;
rho             = obj.HeatTransfer_params.HFS.rho;
cp              = obj.HeatTransfer_params.HFS.cp;
Aboard          = obj.HeatTransfer_params.HFS.A;
epsilonboard    = obj.HeatTransfer_params.HFS.epsilon;

if (strcmp(obj.HeatTransfer_params.HFS.Type,'PCB'))
    kplatex = obj.HeatTransfer_params.HFS.lambdax*s;
    kplatey = obj.HeatTransfer_params.HFS.lambday*s;
else
    kplatex = obj.HeatTransfer_params.HFS.k;
    kplatey = obj.HeatTransfer_params.HFS.k;
end

g               = obj.HeatTransfer_params.constants.g;
sigma           = obj.HeatTransfer_params.constants.sigma;

L_char          = obj.HeatTransfer_params.conditions.L_char;
I               = obj.HeatTransfer_params.conditions.I;
V               = obj.HeatTransfer_params.conditions.V;
Tamb_cold       = obj.HeatTransfer_params.conditions.Tamb(1); 
Tamb_hot        = obj.HeatTransfer_params.conditions.Tamb(2);

threshold = 100; % Rough estimation for minimum temperature in Celsius degrees:
Thot      = checkCelsius(Thot,threshold);
Tcold     = checkCelsius(Tcold,threshold);
Tamb_hot  = checkCelsius(Tamb_hot,threshold);
Tamb_cold = checkCelsius(Tamb_cold,threshold);

errorT       = obj.HeatTransfer_params.Error.errorT;
errorTamb    = obj.HeatTransfer_params.Error.errorTamb;
errorV       = obj.HeatTransfer_params.Error.errorV;
errorI       = obj.HeatTransfer_params.Error.errorI;
errorEpsilon = obj.HeatTransfer_params.Error.errorEpsilon;
errorrho     = obj.HeatTransfer_params.Error.errorrho;
errorcp      = obj.HeatTransfer_params.Error.errorcp;
errors       = obj.HeatTransfer_params.Error.errors;
errorAboard  = obj.HeatTransfer_params.Error.errorAboard;
errorkplate  = obj.HeatTransfer_params.Error.errorkplate;
errorLchar   = obj.HeatTransfer_params.Error.errorLchar;

if obj.HeatTransfer_params.compute_St
    err=fieldnames(obj.HeatTransfer_params.Error);
    if ~any(strcmp(err,'errorUinf'))
        error(['PIRT:parse_PIRT_HeatTransfer: The error in the Uinf ' ...
            'value must be introduced to calculate the heat transfer ' ...
            'error estimation for St'])
    end
    if ~any(strcmp(err,'errorrhoinf'))
        error(['PIRT:parse_PIRT_HeatTransfer: The error in the rhoinf ' ...
            'value must be introduced to calculate the heat transfer ' ...
            'error estimation for St'])
    end
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    errorUinf = obj.HeatTransfer_params.Error.errorUinf;
    rhoinf    = obj.HeatTransfer_params.conditions.rhoinf;
    errorrhoinf = obj.HeatTransfer_params.Error.errorrhoinf;
end

%% Compute heat transfer
fields = fieldnames(obj.result);
if ~any(strcmp(fields,'h'))
    prev = obj.HeatTransfer_params.compute_h;
    obj.HeatTransfer_params.compute_h = 1;
    obj = obj.Calculate_HeatTransfer();
    obj.HeatTransfer_params.compute_h = prev;
    clear prev
    warning(['PIRT:Calculate_HeatTransfer_Error: The h coefficient ' ...
            'must be computed to compute the error by Moffat'])
end
%% Calculating Errors
% Joule losses
sigma_qj = sqrt((errorV.*I/Aboard).^2 + ...
                (errorI.*V./Aboard).^2 + ...
                (errorAboard.*V*I/(Aboard^2)).^2);

% Radiation losses
sigma_qr = sqrt((errorEpsilon.*sigma.*(Thot.^4-Tamb_hot.^4)).^2 + ...
                (errorT.*4.*sigma.*epsilonboard.*Thot.^3).^2 + ...
                (errorTamb.*4.*sigma.*epsilonboard.*Tamb_hot.^3).^2 );

% Tangential losses
if obj.HeatTransfer_params.spatial_der
    fields = fieldnames(obj.result);
    if and(any(strcmp(fields,'d2Tdx2_hot')),any(strcmp(fields,'d2Tdy2_hot')))
        fields = fieldnames(obj.HeatTransfer_params.conditions);
        if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
            fields = fieldnames(obj.HeatTransfer_params.Error);
            if any(strcmp(fields,'errordx'))&&any(strcmp(fields,'errordy'))

                errordx     = obj.HeatTransfer_params.Error.errordx;
                errordy     = obj.HeatTransfer_params.Error.errordy;
                d2Tdx2      = obj.result.d2Tdx2_hot;
                d2Tdy2      = obj.result.d2Tdy2_hot;
                dx          = obj.HeatTransfer_params.conditions.dx;
                dy          = obj.HeatTransfer_params.conditions.dy;

                sigma_d2Tdx2 = sqrt(6*errorT^2/(dx^4) + ...
                    8*errordx^2*(d2Tdx2./dx).^2);
                sigma_d2Tdy2 = sqrt(6*errorT^2/(dy^4) + ...
                    8*errordy^2*(d2Tdy2./dy).^2);
                sigma_qk = sqrt((errorkplate.*d2Tdx2).^2+...
                    (errorkplate.*d2Tdy2).^2+ ...
                    (sigma_d2Tdx2.*kplatex).^2+...
                    (sigma_d2Tdy2.*kplatey).^2);
                clear sigma_d2Tdy2 sigma_d2Tdx2

            else
                sigma_qk = 0;
                warning(['PIRT:Calculate_HeatTransfer_Error: The error in the ' ...
                    'd2Tdx2 and d2Tdy2 can only be computed if the derivatives'...
                    ' were computed with the gradient, errordx and errordy must be '...
                    'included in the conditions.'])
            end
        else
            sigma_qk = 0;
            warning(['PIRT:Calculate_HeatTransfer_Error: The error in the ' ...
                'd2Tdx2 and d2Tdy2 can only be computed if the derivatives'...
                ' were computed with the gradient, dx and dy must be '...
                'included in the conditions.'])
        end

    else
        sigma_qk = 0;
        warning(['PIRT:Calculate_HeatTransfer_Error: To add the spatial derivative ' ...
            'terms to the heat transfer error calculation' ...
            'the temporal derivative has to be computed previously'...
            ' by the gradient method'])
    end
else
    sigma_qk = 0;
end

% Unsteady losses
if obj.HeatTransfer_params.time_der
    fields = fieldnames(obj.result);
    if any(strcmp(fields,'dTdt_hot'))
        fields = fieldnames(obj.HeatTransfer_params.conditions);
        if strcmp(fields,'dt')
            dTdt = obj.result.dTdt_hot;
            dt   = obj.HeatTransfer_params.conditions.dt;

            sigma_qt = sqrt((errorrho*cp*s*dTdt).^2 +...
                            (errorcp*rho*s*dTdt).^2 +...
                            (errors*rho*cp*dTdt).^2 +...
                            (errorT/(sqrt(2)*dt)*rho*cp*s).^2 );

        else
            sigma_qt = 0;
            warning(['PIRT:Calculate_HeatTransfer_Error: The error in the ' ...
                'dTdt can only be computed if the derivatives'...
                ' were computed with the gradient, dt must be '...
                'included in the conditions.'])
        end
    else
        sigma_qt = 0;
        warning(['PIRT:Calculate_HeatTransfer_Error: To add the unsteady ' ...
            'terms to the heat transfer error calculation' ...
            'the temporal derivative has to be computed previously'...
            ' by the gradient method or by applying the sgolay32 filter'])
    end
else
    sigma_qt = 0;
end

% Custom q error
fields = fieldnames(obj.HeatTransfer_params.Error);
if any(strcmp(fields,'sigma_qcustom'))
    sigma_qcustom = obj.HeatTransfer_params.Error.sigma_qcustom;
else
    sigma_qcustom = 0;
end

%% Compute error

sigma_h = sqrt((sigma_qj.^2 + sigma_qr.^2 + sigma_qk.^2 + sigma_qt.^2 + ...
              sigma_qcustom.^2)./(Thot - Tcold*Tamb_hot/Tamb_cold).^2  + ...
              2*errorT^2*obj.result.h.^2);

if obj.HeatTransfer_params.compute_h
    obj.result.errorh = sigma_h;
end

if obj.HeatTransfer_params.compute_Nu
    Tcold_ref = mean(Tcold,3);
    Tfilm     = (Thot + Tcold_ref)/2;
    k         = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; % Air thermal conductivity as a function of T
    
    sigma_k = errorT.*(3*1.5207E-11.*Tfilm.^2 - 2*4.8574E-08.*Tfilm+1.0184E-04)./sqrt(2);
    sigma_Nu = sqrt((sigma_h.*L_char./k).^2 +...
                    (errorLchar.*h./k).^2 +...
                    (sigma_k.*h.*L_char./(k.^2)).^2);
    clear sigma_k
    obj.result.errorNu = sigma_Nu;
end

if obj.HeatTransfer_params.compute_St
    sigma_St = sqrt(sigma_h.^2 + ...
                    errorrhoinf^2./rhoinf.^2+...
                    errorcp.^2./cp.^2 + ...
                    errorUinf.^2./Uinf.^2)/(rhoinf*cp*Uinf);
    obj.result.errorSt = sigma_St;
end

end


function obj = montecarlo_method(obj)
%% Obtain data

if ~isempty(obj.result)
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            if isfile(strcat([obj.output.path,'Thot_filtered.mat']))
                load(strcat([obj.output.path,'Thot_filtered.mat']))
            else
                error('PIRT:Calculate_HeatTransfer: The Thot matrix could not be found')
            end
        else
            if isfile(strcat([obj.output.path,'Thot_filtered.mat']))
                load(strcat([obj.output.path,'\Thot_filtered.mat']))
            else
                error('PIRT:Calculate_HeatTransfer: The Thot matrix could not be found')
            end
        end
        Tcold = obj.Tcold;
    else
        if isfield(obj.result,'Thot_new')
            Thot = obj.result.Thot_new;
        else
            Thot = obj.Thot;
        end
        if isfield(obj.result,'Tcold_new')
            Tcold = obj.result.Tcold_new;
        else
            Tcold = obj.Tcold;
        end
    end
else
    Thot = obj.Thot;
    Tcold = obj.Tcold;
end

Thot = mean(Thot,3);
Tcold = mean(Tcold,3);

s               = obj.HeatTransfer_params.HFS.s;
rho             = obj.HeatTransfer_params.HFS.rho;
cp              = obj.HeatTransfer_params.HFS.cp;
Aboard          = obj.HeatTransfer_params.HFS.A;
epsilonboard    = obj.HeatTransfer_params.HFS.epsilon;

if (strcmp(obj.HeatTransfer_params.HFS.Type,'PCB'))
    kplatex = obj.HeatTransfer_params.HFS.lambdax*s;
    kplatey = obj.HeatTransfer_params.HFS.lambday*s;
else
    kplatex = obj.HeatTransfer_params.HFS.k;
    kplatey = obj.HeatTransfer_params.HFS.k;
end

g               = obj.HeatTransfer_params.constants.g;
sigma           = obj.HeatTransfer_params.constants.sigma;

L_char          = obj.HeatTransfer_params.conditions.L_char;
I               = obj.HeatTransfer_params.conditions.I;
V               = obj.HeatTransfer_params.conditions.V;
Tamb_cold       = obj.HeatTransfer_params.conditions.Tamb(1); 
Tamb_hot        = obj.HeatTransfer_params.conditions.Tamb(2);

threshold = 100; % Rough estimation for minimum temperature in Celsius degrees:
Thot      = checkCelsius(Thot,threshold);
Tcold     = checkCelsius(Tcold,threshold);
Tamb_hot  = checkCelsius(Tamb_hot,threshold);
Tamb_cold = checkCelsius(Tamb_cold,threshold);

errorT       = obj.HeatTransfer_params.Error.errorT;
errorTamb    = obj.HeatTransfer_params.Error.errorTamb;
errorV       = obj.HeatTransfer_params.Error.errorV;
errorI       = obj.HeatTransfer_params.Error.errorI;
errorEpsilon = obj.HeatTransfer_params.Error.errorEpsilon;
errorrho     = obj.HeatTransfer_params.Error.errorrho;
errorcp      = obj.HeatTransfer_params.Error.errorcp;
errors       = obj.HeatTransfer_params.Error.errors;
errorAboard  = obj.HeatTransfer_params.Error.errorAboard;
errorkplate  = obj.HeatTransfer_params.Error.errorkplate;
errorLchar   = obj.HeatTransfer_params.Error.errorLchar;
errork       = obj.HeatTransfer_params.Error.errork;

if obj.HeatTransfer_params.compute_St
    err=fieldnames(obj.HeatTransfer_params.Error);
    if ~any(strcmp(err,'errorUinf'))
        error(['PIRT:parse_PIRT_HeatTransfer: The error in the Uinf ' ...
            'value must be introduced to calculate the heat transfer ' ...
            'error estimation for St'])
    end
    if ~any(strcmp(err,'errorrhoinf'))
        error(['PIRT:parse_PIRT_HeatTransfer: The error in the rhoinf ' ...
            'value must be introduced to calculate the heat transfer ' ...
            'error estimation for St'])
    end
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    errorUinf = obj.HeatTransfer_params.Error.errorUinf;
    rhoinf    = obj.HeatTransfer_params.conditions.rhoinf;
    errorrhoinf = obj.HeatTransfer_params.Error.errorrhoinf;
end

Tfilm = (Thot + Tcold)/2;
k     = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; % Air thermal conductivity as a function of T

err=fieldnames(obj.HeatTransfer_params.Error);
if ~any(strcmp(err,'samples'))
    warning('PIRT:parse_PIRT_HeatTransfer: Default number of samples taken as 1000. The number of samples was not introduced for the Montecarlo error estimation')
    n = 1000;
else
    n = obj.HeatTransfer_params.Error.samples;
end

tder = 1;
if strcmp(obj.output.type,'file')
    if endsWith(obj.output.path,'\')
        if isfile(strcat([obj.output.path,'dTdt.mat']))
            load(strcat([obj.output.path,'dTdt.mat']))
            tder = 1;
        else
            tder = 0;
        end
    else
        if isfile(strcat([obj.output.path,'\dTdt.mat']))
            load(strcat([obj.output.path,'\dTdt.mat']))
            tder = 1;
        else
            tder = 0;
        end
    end
    dTdt_m = mean(dTdt_hot,3);
    clear dTdt_hot
else
    fields = fieldnames(obj.result);
    if any(strcmp(fields,'dTdt_hot'))
        dTdt_m = mean(obj.result.dTdt_hot,3);
        tder = 1;
    else
        tder = 0;
    end
end

sder = 1;
if strcmp(obj.output.type,'file')
    if endsWith(obj.output.path,'\')
        if isfile(strcat([obj.output.path,'d2Tdx2.mat']))&&isfile(strcat([obj.output.path,'d2Tdy2.mat']))
            load(strcat([obj.output.path,'d2Tdx2.mat']))
            load(strcat([obj.output.path,'d2Tdy2.mat']))
            sder = 1;
        else
            sder = 0;
        end
    else
        if isfile(strcat([obj.output.path,'\d2Tdx2.mat']))&&isfile(strcat([obj.output.path,'\d2Tdy2.mat']))
            load(strcat([obj.output.path,'\d2Tdx2.mat']))
            load(strcat([obj.output.path,'\d2Tdy2.mat']))
            sder = 1;
        else
            sder = 0;
        end
    end
    d2Tdx2_m = mean(d2Tdx2_hot,3);
    d2Tdy2_m = mean(d2Tdy2_hot,3);
    clear d2Tdx2_hot d2Tdy2_hot
else
    fields = fieldnames(obj.result);
    if any(strcmp(fields,'d2Tdx2_hot'))&&any(strcmp(fields,'d2Tdy2_hot'))
        d2Tdx2_m = mean(obj.result.d2Tdx2_hot,3);
        d2Tdy2_m = mean(obj.result.d2Tdy2_hot,3);
        tder = 1;
    else
        tder = 0;
    end
end


for i=1:n

    current.Thot      = normrnd(Thot,Thot.*errorT);
    current.Tcold     = normrnd(Tcold,Tcold.*errorT);
    current.Tamb_cold = normrnd(Tamb_cold,Tamb_cold.*errorTamb);
    current.Tamb_hot  = normrnd(Tamb_hot,Tamb_hot.*errorTamb);
    current.V         = normrnd(V,V.*errorV);
    current.I         = normrnd(I,I.*errorI);
    current.epsilon   = normrnd(epsilonboard,epsilonboard.*errorEpsilon);
    current.rho       = normrnd(rho,rho.*errorrho);
    current.cp        = normrnd(cp,cp.*errorcp);
    current.s         = normrnd(s,s.*errors);
    current.kplatex   = normrnd(kplatex,kplatex.*errorkplate);
    current.kplatey   = normrnd(kplatey,kplatey.*errorkplate);
    current.Aboard    = normrnd(Aboard,Aboard.*errorAboard);
    current.L_char    = normrnd(L_char,L_char.*errorLchar);
    current.k         = normrnd(k,k.*errork);

    ratio   = current.Tamb_hot/current.Tamb_cold;

    qj      = current.V*current.I/(current.Aboard); %[W/m^2]

    %-- Radiative heat flux: qr = σ·ε·(Th⁴-Tc⁴)
    qrad    = sigma.*current.epsilon.*(current.Thot.^4-current.Tamb_hot.^4); % [W/m^2]

    %-- Internal energy
    if obj.HeatTransfer_params.time_der&&tder==1
        errordTdt = obj.HeatTransfer_params.Error.errordTdt;
        dTdt_hot  = normrnd(dTdt_m,dTdt_m.*errordTdt);
        unsteady  = current.rho*current.s*current.cp*dTdt_hot; % [W/m^2]
    else
        unsteady=0;
    end

    %-- Tangencial conduction heat flux
    if obj.HeatTransfer_params.spatial_der&&sder==1

        errord2Tdx2 = obj.HeatTransfer_params.Error.errord2Tdx2;
        errord2Tdy2 = obj.HeatTransfer_params.Error.errord2Tdy2;
        d2dTdx2     = normrnd(d2Tdx2_m,d2Tdx2_m.*errord2Tdx2);
        d2dTdy2     = normrnd(d2Tdy2_m,d2Tdy2_m.*errord2Tdy2);
        qk          = current.kplatex.*d2dTdx2+current.kplatey.*d2dTdy2; % [W/m^2]

    else
        qk=0;
    end

    h(:,:,i)  = ( -qk - qrad + qj - unsteady ) ./ (current.Thot - (current.Tcold*ratio) ); %[W/K m^2]
    if obj.HeatTransfer_params.compute_Nu
        Nu(:,:,i) = h(:,:,i).*current.L_char./current.k;
    end
    if obj.HeatTransfer_params.compute_St
        current.Uinf    = normrnd(Uinf,Uinf.*errorUinf);
        current.rhoinf    = normrnd(rhoinf,rhoinf.*errorrhoinf);
        St(:,:,i) = h(:,:,i)./(current.rho*current.cp*current.Uinf);
    end
end

if obj.HeatTransfer_params.compute_h
    errorh   = squeeze(mean(mean(h,1),2));
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'errorh.mat']),"errorh",'-v7.3')
            if isfile(strcat([obj.output.path,'h.mat']))
                hres = load(strcat([obj.output.path,'h.mat']));
                errorh_p = abs(squeeze(mean(mean(h,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'errorh_p.mat']),"errorh_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        else
            save(strcat([obj.output.path,'\errorh.mat']),"errorh",'-v7.3')
            if isfile(strcat([obj.output.path,'\h.mat']))
                hres = load(strcat([obj.output.path,'\h.mat']));
                errorh_p = abs(squeeze(mean(mean(h,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'\errorh_p.mat']),"errorh_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        end
    else
        hres = obj.result.h;
        errorh_p = abs(squeeze(mean(mean(h,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
        obj.result.errorh = errorh;
        obj.result.errorh_p = errorh_p;
    end
    fprintf('The mean value of h is %.2f  with a standard deviation %.2f \n',nanmean(errorh),nanstd(errorh));
end


if obj.HeatTransfer_params.compute_Nu
    errorNu   = squeeze(mean(mean(Nu,1),2));
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'errorNu.mat']),"errorNu",'-v7.3')
            if isfile(strcat([obj.output.path,'Nu.mat']))
                hres = load(strcat([obj.output.path,'Nu.mat']));
                errorNu_p = abs(squeeze(mean(mean(Nu,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'errorNu_p.mat']),"errorNu_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        else
            save(strcat([obj.output.path,'\errorNu.mat']),"errorNu",'-v7.3')
            if isfile(strcat([obj.output.path,'\Nu.mat']))
                hres = load(strcat([obj.output.path,'\Nu.mat']));
                errorNu_p = abs(squeeze(mean(mean(Nu,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'\errorh_p.mat']),"errorNu_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        end
    else
        hres = obj.result.Nu;
        errorNu_p = abs(squeeze(mean(mean(Nu,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
        obj.result.errorNu = errorNu;
        obj.result.errorNu_p = errorNu_p;
    end
    fprintf('The mean value of Nu is %.2f  with a standard deviation %.2f \n',nanmean(errorNu),nanstd(errorNu));
end


if obj.HeatTransfer_params.compute_St
    errorSt   = squeeze(mean(mean(St,1),2));
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'errorSt.mat']),"errorSt",'-v7.3')
            if isfile(strcat([obj.output.path,'St.mat']))
                hres = load(strcat([obj.output.path,'St.mat']));
                errorSt_p = abs(squeeze(mean(mean(St,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'errorSt_p.mat']),"errorSt_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        else
            save(strcat([obj.output.path,'\errorSt.mat']),"errorSt",'-v7.3')
            if isfile(strcat([obj.output.path,'\St.mat']))
                hres = load(strcat([obj.output.path,'\St.mat']));
                errorSt_p = abs(squeeze(mean(mean(St,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
                save(strcat([obj.output.path,'\errorh_p.mat']),"errorSt_p",'-v7.3')
            else
                warning('PIRT:Calculate_HeatTransfer_Error: The percentage error could not be computed since the h data was not found.')
            end
        end
    else
        hres = obj.result.St;
        errorSt_p = abs(squeeze(mean(mean(St,1),2))-mean(hres,'all'))./mean(hres,'all').*100;
        obj.result.errorSt = errorSt;
        obj.result.errorSt_p = errorSt_p;
    end
    fprintf('The mean value of St is %.2f  with a standard deviation %.2f \n',nanmean(errorSt),nanstd(errorSt));
end

end

% Aux functions:
function T = checkCelsius(T,threshold)
cel2kel = 273.15;
Tm = mean(T,'all','omitnan');
if Tm < threshold
    T = T + cel2kel;
    warning('PIRT:CalculateHeatTransfer: Temperature was converted to Kelvin. Check your inputs!')
end
end