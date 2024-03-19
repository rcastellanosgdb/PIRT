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
kplate          = obj.HeatTransfer_params.HFS.k;
Aboard          = obj.HeatTransfer_params.HFS.A;
epsilonboard    = obj.HeatTransfer_params.HFS.epsilon;

g               = obj.HeatTransfer_params.constants.g;
sigma           = obj.HeatTransfer_params.constants.sigma;

L_char          = obj.HeatTransfer_params.conditions.L_char;
I               = obj.HeatTransfer_params.conditions.I;
V               = obj.HeatTransfer_params.conditions.V;
Tamb_cold       = obj.HeatTransfer_params.conditions.Tamb(1) + 273.15; % Convert to Kelvin
Tamb_hot        = obj.HeatTransfer_params.conditions.Tamb(2) + 273.15;

errorTaw = obj.HeatTransfer_params.Error.errorT;
errorTw = obj.HeatTransfer_params.Error.errorT;
errorTamb = obj.HeatTransfer_params.Error.errorTamb;
errorV = obj.HeatTransfer_params.Error.errorV;
errorI = obj.HeatTransfer_params.Error.errorI;
errorEpsilon = obj.HeatTransfer_params.Error.errorEpsilon;
errorrho = obj.HeatTransfer_params.Error.errorrho;
errorcp = obj.HeatTransfer_params.Error.errorcp;
errors = obj.HeatTransfer_params.Error.errors;
errorAboard = obj.HeatTransfer_params.Error.errorAboard;
errorkplate = obj.HeatTransfer_params.Error.errorkplate;
errorLchar = obj.HeatTransfer_params.Error.errorLchar;
errork = obj.HeatTransfer_params.Error.errork;

if obj.HeatTransfer_params.compute_St
    err=fieldnames(obj.HeatTransfer_params.Error);
    if ~any(strcmp(err,'errorUinf'))
        error(['PIRT:parse_PIRT_HeatTransfer: The error in the Uinf ' ...
            'value must be introduced to calculate the heat transfer ' ...
            'error estimation for St'])
    end
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    errorUinf = obj.HeatTransfer_params.Error.errorUinf;
end


Tcold_ref = mean(Tcold,3);
Tfilm   = (Thot + Tcold_ref)/2;
k = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; % Air thermal conductivity as a function of T
%% Calculating derivatives
%Joule losses
dNudQj    = L_char./(k.*(Thot-Tcold));
dQjdV     = I/Aboard;
dQjdI     = V/Aboard;
dQjdAs     = -1*V*I/Aboard^2;
%Radiation losses
dNudQr    = -L_char./(k.*(Thot-Tcold));
dQrdeps   = sigma.*(Thot.^4-Tamb_hot.^4);
dQrdTw    = 4.*sigma.*epsilonboard.*Thot.^3;
dQrdTamb  = -4.*sigma.*epsilonboard.*Tamb_hot.^3;

%% Calculating errors

errorQj   = ( (dQjdV.*errorV).^2 + (dQjdI.*errorI).^2 + (dQjdAs.*errorAboard).^2 ).^0.5;
errorQr   = ( (dQrdeps.*errorEpsilon).^2 + (dQrdTw.*errorTw).^2 + (dQrdTamb.*errorTamb).^2 ).^0.5;

if obj.HeatTransfer_params.time_der
    fields = fieldnames(obj.result);
    if any(strcmp(fields,'dTdt_hot'))
        fields = fieldnames(obj.HeatTransfer_params.Error);
        if any(strcmp(fields,'errordTdt'))
            errordTdt = obj.HeatTransfer_params.Error.errordTdt;

            dTdt = obj.result.dTdt_hot;

            dNudQuns  = -L_char./(k.*(Thot-Tcold));
            dQunsdcp  = rho.*s.*dTdt;
            dQunsdrho = cp.*s.*dTdt;
            dQunsds   = cp.*rho.*dTdt;
            dQunsdTs  = cp.*rho.*s;

            errorQuns = ( (dQunsdcp.*errorcp).^2 + (dQunsdrho.*errorrho).^2 + (dQunsds.*errors).^2 + (dQunsdTs.*errordTdt).^2 ).^0.5;
        else
            dNudQuns = 0;
            errorQuns = 0;
            warning(['PIRT:Calculate_HeatTransfer_Error: The error in the ' ...
                'dTdt value must be introduced to calculate the heat ' ...
                'transfer error estimation with the time derivative option'])

        end
    else
        dNudQuns = 0;
        errorQuns = 0;
        warning(['PIRT:Calculate_HeatTransfer_Error: To add the unsteady ' ...
            'terms to the heat transfer error calculation an sgolay32 ' ...
            'filter has to be performed on the hot images'])
    end
else
    errorQuns = 0;
end

if obj.HeatTransfer_params.spatial_der
    fields = fieldnames(obj.result);
    if and(any(strcmp(fields,'d2Tdx2_hot')),any(strcmp(fields,'d2Tdy2_hot')))
        fields = fieldnames(obj.HeatTransfer_params.Error);
        if and(any(strcmp(fields,'errord2Tdx2')),any(strcmp(fields,'errord2Tdx2')))
            errord2Tdx2 = obj.HeatTransfer_params.Error.errord2Tdx2;
            errord2Tdy2 = obj.HeatTransfer_params.Error.errord2Tdy2;
            d2Tdx2 = obj.result.d2Tdx2_hot;
            d2Tdy2 = obj.result.d2Tdy2_hot;

            dNudQk    = -L_char./(k.*(Thot-Tcold));
            dQkdkcirc = s.*(d2Tdx2 + d2Tdy2);
            dQkds     = kplate.*(d2Tdx2 + d2Tdy2);
            dQdT2     = kplate.*s;

            errorQk   = ( (dQkdkcirc.*errorkplate).^2 + (dQkds*errors).^2 + (dQdT2.*errord2Tdx2).^2 + (dQdT2.*errord2Tdy2).^2).^0.5;
        else
            dNudQk = 0;
            errorQk = 0;
            error(['PIRT:Calculate_HeatTransfer_Error: The error in the ' ...
                'd2dTdx2 and d2Tdy2 values must be introduced to calculate ' ...
                'the heat transfer error estimation with the spatial ' ...
                'derivative option'])
        end
    else
        dNudQk = 0;
        errorQk = 0;
        warning(['PIRT:Calculate_HeatTransfer_Error: To add the spatial ' ...
            'derivative terms to the heat transfer error calculation ' ...
            'an sgolay32 filter has to be performed on the hot images'])
    end
else
    dNudQk = 0;
    errorQk = 0;
end

%- Independent terms:
fields = fieldnames(obj.result);
if any(strcmp(fields,'Nu'))
    Nu = obj.result.Nu;
    h = Nu.*k./L_char;
    dNudk     = -Nu./k;
    dNudLchar = +Nu./L_char;
    dNudTw    = -Nu./(Thot-Tcold);
    dNudTaw   = +Nu./(Thot-Tcold);
else
    error('PIRT:Calculate_HeatTransfer_Error: The Nusselt number must be computed before trying to estimate the error')
end


errorh = ((dNudQj.*errorQj).^2 + (dNudQr.*errorQr).^2 + (dNudQk.*errorQk).^2 + (dNudQuns.*errorQuns).^2 +...
    (dNudTw.*errorTw).^2 + (dNudTaw.*errorTaw).^2 ).^0.5;
errorh_p = errorh./h.*100;

errorNu   = errorh + (dNudk.*errork).^2 + (dNudLchar.*errorLchar).^2;
errorNu_p = errorNu./Nu.*100;

if obj.HeatTransfer_params.compute_h
    obj.result.errorh      = errorh;
    obj.result.errorh_p = errorh_p;
end

if obj.HeatTransfer_params.compute_Nu
    obj.result.errorNu      = errorNu;
    obj.result.errorNu_p = errorNu_p;
end

if obj.HeatTransfer_params.compute_St
    dStdrho = h./(cp*Uinf);
    dStdcp= h./(rho*Uinf);
    dStdUinf = h./(rho*cp);
    obj.result.errorSt = errorh + (dStdrho*errorrho).^2 + (dStdcp*errorcp).^2 + (dStdUinf*errorUinf).^2;
    St = h./(rho*cp*Uinf);
    obj.result.errorSt_p = obj.result.errorSt./St.*100;
end

end


function obj = montecarlo_method(obj)
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

Thot = mean(Thot,3);
Tcold = mean(Tcold,3);

s               = obj.HeatTransfer_params.HFS.s;
rho             = obj.HeatTransfer_params.HFS.rho;
cp              = obj.HeatTransfer_params.HFS.cp;
kplate          = obj.HeatTransfer_params.HFS.k;
Aboard          = obj.HeatTransfer_params.HFS.A;
epsilonboard    = obj.HeatTransfer_params.HFS.epsilon;

g               = obj.HeatTransfer_params.constants.g;
sigma           = obj.HeatTransfer_params.constants.sigma;

L_char          = obj.HeatTransfer_params.conditions.L_char;
I               = obj.HeatTransfer_params.conditions.I;
V               = obj.HeatTransfer_params.conditions.V;
Tamb_cold       = obj.HeatTransfer_params.conditions.Tamb(1) + 273.15; % Convert to Kelvin
Tamb_hot        = obj.HeatTransfer_params.conditions.Tamb(2) + 273.15;

errorTaw = obj.HeatTransfer_params.Error.errorT;
errorTw = obj.HeatTransfer_params.Error.errorT;
errorTamb = obj.HeatTransfer_params.Error.errorTamb;
errorV = obj.HeatTransfer_params.Error.errorV;
errorI = obj.HeatTransfer_params.Error.errorI;
errorEpsilon = obj.HeatTransfer_params.Error.errorEpsilon;
errorrho = obj.HeatTransfer_params.Error.errorrho;
errorcp = obj.HeatTransfer_params.Error.errorcp;
errors = obj.HeatTransfer_params.Error.errors;
errorAboard = obj.HeatTransfer_params.Error.errorAboard;
errorkplate = obj.HeatTransfer_params.Error.errorkplate;
errorLchar = obj.HeatTransfer_params.Error.errorLchar;
errork = obj.HeatTransfer_params.Error.errork;
if obj.HeatTransfer_params.compute_St
    err=fieldnames(obj.HeatTransfer_params.Error);
    if ~any(strcmp(err,'errorUinf'))
        error('PIRT:parse_PIRT_HeatTransfer: The error in the Uinf value must be introduced to calculate the heat transfer error estimation for St')
    end
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    errorUinf = obj.HeatTransfer_params.Error.errorUinf;
end


Tfilm   = (Thot + Tcold)/2;
k = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; % Air thermal conductivity as a function of T

err=fieldnames(obj.HeatTransfer_params.Error);
if ~any(strcmp(err,'samples'))
    warning('PIRT:parse_PIRT_HeatTransfer: Default number of samples taken as 1000. The number of samples was not introduced for the Montecarlo error estimation')
    n=1000;
else
    n = obj.HeatTransfer_params.Error.samples;
end

for i=1:n

    current.Thot = normrnd(Thot,Thot.*errorTw);
    current.Tcold = normrnd(Tcold,Tcold.*errorTaw);
    current.Tamb_cold = normrnd(Tamb_cold,Tamb_cold.*errorTamb);
    current.Tamb_hot = normrnd(Tamb_hot,Tamb_hot.*errorTamb);
    current.V = normrnd(V,V.*errorV);
    current.I = normrnd(I,I.*errorI);
    current.epsilon = normrnd(epsilonboard,epsilonboard.*errorEpsilon);
    current.rho = normrnd(rho,rho.*errorrho);
    current.cp = normrnd(cp,cp.*errorcp);
    current.s = normrnd(s,s.*errors);
    current.kplate = normrnd(kplate,kplate.*errorkplate);
    current.Aboard = normrnd(Aboard,Aboard.*errorAboard);
    current.L_char = normrnd(L_char,L_char.*errorLchar);
    current.k = normrnd(k,k.*errork);

    ratio   = current.Tamb_hot/current.Tamb_cold;

    qj      = current.V*current.I/(current.Aboard); %[W/m^2]

    %-- Radiative heat flux: qr = σ·ε·(Th⁴-Tc⁴)
    qrad    = sigma.*current.epsilon.*(current.Thot.^4-current.Tamb_hot.^4); % [W/m^2]

    %-- Internal energy
    if obj.HeatTransfer_params.time_der
        if isempty(obj.filter_params)
            unsteady = 0;
        elseif any(cellfun(@isempty,obj.filter_params.sgolay32_params))
            fields = fieldnames(obj.result);
            if any(strcmp(fields,'dTdt_hot'))
                errordTdt = obj.HeatTransfer_params.Error.errordTdt;
                dTdt_hot = normrnd(mean(obj.result.dTdt_hot,3),mean(obj.result.dTdt_hot,3).*errordTdt);
                unsteady = current.rho*current.s*current.cp*dTdt_hot; % [W/m^2]
            else
                unsteady = 0;
            end
        else
            unsteady = 0;
        end
    else
        unsteady=0;
    end

    %-- Tangencial conduction heat flux
    if obj.HeatTransfer_params.spatial_der
        if isempty(obj.filter_params)
            qk=0;
        elseif any(cellfun(@isempty,obj.filter_params.sgolay32_params))
            fields = fieldnames(obj.result);
            if and(any(strcmp(fields,'d2Tdx2_hot')),any(strcmp(fields,'d2Tdy2_hot')))
                errord2Tdx2 = obj.HeatTransfer_params.Error.errord2Tdx2;
                errord2Tdy2 = obj.HeatTransfer_params.Error.errord2Tdy2;
                d2dTdx2 = normrnd(mean(obj.result.d2Tdx2_hot,3),mean(obj.result.d2Tdx2_hot,3).*errord2Tdx2);
                d2dTdy2 = normrnd(mean(obj.result.d2Tdy2_hot,3),mean(obj.result.d2Tdy2_hot,3).*errord2Tdy2);
                qk = current.kplate*current.s*(d2dTdx2+d2dTdy2); % [W/m^2]
            else
                qk = 0;
            end
        else
            qk = 0;
        end
    else
        qk=0;
    end

    h(:,:,i)  = ( -qk - qrad + qj - unsteady ) ./ (current.Thot - (current.Tcold*ratio) ); %[W/K m^2]
    if obj.HeatTransfer_params.compute_Nu
        Nu(:,:,i) = h(:,:,i).*current.L_char./current.k;
    end
    if obj.HeatTransfer_params.compute_St
        current.Uinf    = normrnd(Uinf,Uinf.*errorUinf);
        St(:,:,i) = h(:,:,i)./(current.rho*current.cp*current.Uinf);
    end
end

if obj.HeatTransfer_params.compute_h
    obj.result.errorh      = mean(h,3);
    obj.result.errorh_p = mean(h,3)./mean(obj.result.h,3).*100;
end

if obj.HeatTransfer_params.compute_Nu
    obj.result.errorNu      = mean(Nu,3);
    obj.result.errorNu_p = mean(Nu,3)./mean(obj.result.Nu,3).*100;
end

if obj.HeatTransfer_params.compute_St
    obj.result.errorSt      = mean(St,3);
    obj.result.errorSt_p = mean(St,3)./mean(obj.result.St,3).*100;
end

end