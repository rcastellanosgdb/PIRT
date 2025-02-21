function obj = Calculate_HeatTransfer(obj)
%CALCULATE_HEATTRANSFER is a method of the PIRT class utilized to calculate the
%   nusselt number,Stanton number and convection coefficient. There are no
%   specific inputs since attributes previusly
%   set are utilized. The convection coefficient is computed via an energy
%   balance h = (qj" - qr" - qk" - qb")/(Tw - Taw) where:
%       qj": is the heat flux provided to the heat flux sensor by joule
%           effect
%       qr": radiation heat flux
%       qk": tangential heat flux throufh the foil sensor
%       qr": heat flux through the thin layer of air behind the sensor
%       Tw: wall temperature. Considered as the Thot previously seted (if
%       it is a 3D matrix the third dimension is assumed to be the temporal
%       dimension and it is averaged
%       Taw: adiabatic wall temperature. Which is the Tcold matrix
%       (averaged over the third dimension if there is any) multiplied by
%       the ratio of hot to cold ambient temperatures.
%   The required attributes of the object are:
%       Thot: matrix of temperatrues where the sensor foil was heated
%           through the joule effect
%       Tcold: matrix of temperatures whith no heating to the sensor
%       HFS: struct containing the heat foil sensor parameters:
%           s: foil sensor thickness [m]
%           cp: heat capacity of foil [J/kgK]
%           k: thermal conductivity foil [m K]
%           rho: foil sensor density [kg/m3]
%           Aboard: Area of foil sensor [m2]
%           epsilonboard: foil sensor emissivity [-]
%           sigma: Stefan-Boltzmann constant [W/(m^2 K^4)]
%           g: gravity acceleration [m/s2]
%           Type: Foil or PCB model. Foil is selected as default.
%           lambdax: thermal conductivity in the x direction. Only required
%               in PCB model[m K]
%           lambday: thermal conductivity in the y direction. Only required
%               in PCB model[m K]
%       conditions: struct containing information about the measuring
%       environment
%           L_char: foil sensor characteristic length [m]
%           T_amb: vector of 2 elements containing first the cold ambient
%           temperature and then the hot [ºC or K]
%           V: voltage supplied to the foil sensor [V]
%           I: current supplied to foil sensor [A]
%   A more advanced option is the opportunity to take into account the
%   spatial variations of the temperature and the unsteady terms. However a
%   previosuly computed sgolay32 filter is required to compute the
%   temperature derivatives.
%   Returns:
%       obj: the object with the new filtered images in the obj.result
%       attribute

%   References:
%     [1] T. Astarita, G. M. Carlomagno, Infrared thermography for thermo-
%         fluid-dynamics, Springer Science & Business Media, 2012

%   Author(s): I. Robledo, R. Castellanos
%   Copyright 2023 Universidad Carlos III de Madrid

[Thot,Tcold,s,s_p,rho,rho_p,cp,cp_p,k,kp,A,epsilon,L_char,sigma,PCB,sides,customQ] = parseinputs(obj);

disp('-- Calculating heat transfer maps');
% Check temperatures (convert to Kelvin if the average temperature is lower
% than an expected threshold):
threshold   = 100; % Rough estimation for minimum temperature in Celsius degrees:
Thot        = checkCelsius(Thot,threshold);
Tcold       = checkCelsius(Tcold,threshold);
Tamb_cold   = checkCelsius(obj.HeatTransfer_params.conditions.Tamb(1),threshold);
Tamb_hot    = checkCelsius(obj.HeatTransfer_params.conditions.Tamb(2),threshold);

ratio   = Tamb_hot/Tamb_cold; %[-]

% Initialize heat value
q = 0;

%-- Heat Flux due to Joule Effect: qj = VI/A
q      = q + obj.HeatTransfer_params.conditions.V*obj.HeatTransfer_params.conditions.I/(A); %[W/m^2]

%-- Radiative heat flux: qr = σ·ε·(Th⁴-Tc⁴)
q      = q - sides*sigma.*epsilon.*(Thot.^4-Tamb_hot.^4); % [W/m^2]

%-- Unsteady terms
if obj.HeatTransfer_params.time_der
    [obj,q] = compute_unsteady(q,Thot,obj,rho,rho_p,s,s_p,cp,cp_p);

end

%-- Tangencial conduction heat flux
if obj.HeatTransfer_params.spatial_der
    [obj,q] = compute_tangential(q,Thot,obj,s,k,kp,PCB);
end

if ~isempty(customQ)
    for i=1:length(customQ)
        checkSizeMats(Thot,customQ{i});
        q=q+customQ{i};
    end
end

% Calculate h (convective heat transfer coefficient) h [W/Km2]
h      = q./ ( Thot - (Tcold*ratio) ); %[W/K m^2]

if obj.HeatTransfer_params.compute_h
    if strcmp(obj.output.type,'file')
        disp('--> Saving h information into file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'h.mat']),"h",'-v7.3')
        else
            save(strcat([obj.output.path,'\h.mat']),"h",'-v7.3')
        end
    else
        obj.result.h = h;
    end
end

if obj.HeatTransfer_params.compute_Nu
    % Calculation of Heated Thin-Foil Sensor heat terms:
    %-- Tcold or T adiabatic wall/ T bulk temperature
    %   Compute the film temperature. (Check if available memory is enough
    %   to compute it directly; otherwise, make it in different blocks
    %   iteratively.
    try
        [~,sys]     = memory;
        byte2bit    = 8;
        max_Nimg    = floor(sys.PhysicalMemory.Available /(2*numel(Tcold_ref)*byte2bit));
        Tfilm   = zeros(size(Tcold_ref));
        for n = 0:max_Nimg:size(Thot,3)
            if (n + max_Nimg < size(Thot,3))
                nrange = n + (1:max_Nimg);
            else
                nrange = (n+1):size(Thot,3);
            end
            Tfilm(:,:,nrange)   = (Thot(:,:,nrange)+Tamb_hot)/2; %[K]
        end
    catch
        Tfilm = (Thot+Tamb_hot)/2;
    end
    % Air thermal conductivity as a function of T
    kair   = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; clear Tfilm
    Nu = h.*L_char./kair;

    if strcmp(obj.output.type,'file')
        disp('--> Saving Nu information into file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'Nu.mat']),"Nu",'-v7.3')
        else
            save(strcat([obj.output.path,'\Nu.mat']),"Nu",'-v7.3')
        end
    else
        obj.result.Nu = Nu;
    end
    clear Nu
end

if obj.HeatTransfer_params.compute_St
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    rhoinf  = obj.HeatTransfer_params.conditions.rhoinf;
    St      = h./(rhoinf*cp*Uinf);

    if strcmp(obj.output.type,'file')
        disp('--> Saving St information into file')
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'St.mat']),"St",'-v7.3')
        else
            save(strcat([obj.output.path,'\St.mat']),"St",'-v7.3')
        end
    else
        obj.result.St = St;
    end
    clear St
end
end

%-------
% Aux functions:
function T = checkCelsius(T,threshold)
cel2kel = 273.15;
Tm = mean(T,'all','omitnan');
if Tm < threshold
    T = T + cel2kel;
    warning('PIRT:CalculateHeatTransfer: Temperature was converted to Kelvin. Check your inputs!')
end
end

function [Thot,Tcold,s,s_paint,rho,rho_paint,cp,cp_paint,k,kp,A,epsilon,L_char,sigma,PCB,sides,customQ] = parseinputs(obj)

HFS         = obj.HeatTransfer_params.HFS;
Conditions  = obj.HeatTransfer_params.conditions;
params      = obj.HeatTransfer_params;

if ~isempty(obj.result)
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            if isfile(strcat([obj.output.path,'Thot_filtered.mat']))
                load(strcat([obj.output.path,'Thot_filtered.mat']))
            else
                error('PIRT:Calculate_HeatTransfer: The Thot matrix could not be found')
            end
        else
            if isfile(strcat([obj.output.path,'\Thot_filtered.mat']))
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

if (length(size(Thot))==3)&&(length(size(Tcold))==3)
    if sum(size(Thot)-size(Tcold))~=0
        error('PIRT:Calculate_HeatTransfer: Incompatible images size')
    end
else
    if sum(size(Thot,[1,2])-size(Tcold,[1,2]))~=0
        error('PIRT:Calculate_HeatTransfer: Incompatible images size')
    end
end

if or(isempty(HFS),isempty(Conditions))
    error('PIRT:Calculate_HeatTransfer: Both the HFS and Conditions data must be introduced to calculate the heat transfer')
end

hfs_data = fieldnames(HFS);

if ~any(strcmp(hfs_data,'s'))
    error('PIRT:Calculate_HeatTransfer: the TFS thickness s must be introduced to compute the heat transfer')
end
if ~any(strcmp(hfs_data,'rho'))
    error('PIRT:Calculate_HeatTransfer: rho must be introduced to compute the heat transfer')
end
if ~any(strcmp(hfs_data,'cp'))
    error('PIRT:Calculate_HeatTransfer: cp must be introduced to compute the heat transfer')
end
if ~any(strcmp(hfs_data,'epsilon'))
    error('PIRT:Calculate_HeatTransfer: epsilon must be introduced to compute the heat transfer')
end
if ~any(strcmp(hfs_data,'A'))
    error('PIRT:Calculate_HeatTransfer: either the Area or H and W must be introduced to compute the heat transfer')
end
if ~any(strcmp(hfs_data,'sides'))
    warning('PIRT:Calculate_HeatTransfer: number of exposed sides was not selected, it will be set to a standard value of sides=1')
    sides=1;
else
    sides = obj.HeatTransfer_params.HFS.sides;
end
if ~any(strcmp(hfs_data,'s_paint'))
    warning('PIRT:Calculate_HeatTransfer: the paint thickness s_paint was not introduced,  it will be set to a standard value of s_paint=21.81*10-6 m. Ref. Stafford, Jason / Walsh, Ed / Egan, Vanessa Characterizing convective heat transfer using infrared thermography and the heated-thin-foil technique 2009-09')
    s_paint = 21.81*(10^(-6));
else
    s_paint = obj.HeatTransfer_params.HFS.s_paint;
end
if ~any(strcmp(hfs_data,'rho_paint'))
    warning('PIRT:Calculate_HeatTransfer: the paint density rho_paint was not introduced,  it will be set to a standard value of rho_paint=1300 kg/m3.')
    rho_paint = 1300;
else
    rho_paint = obj.HeatTransfer_params.HFS.rho_paint;
end
if ~any(strcmp(hfs_data,'cp_paint'))
    warning('PIRT:Calculate_HeatTransfer: the paint heat capacity cp_paint was not introduced,  it will be set to a standard value of s_paint=5000 J/kgK.')
    cp_paint = 5000;
else
    cp_paint = obj.HeatTransfer_params.HFS.cp_paint;
end
if ~any(strcmp(hfs_data,'lambda_paint'))
    warning('PIRT:Calculate_HeatTransfer: ther thermal conductivity lambda_paint of the paint was not introduced, it will be set to a standard value of matte black paint lambda_foil=1.38 W/mK. Ref. Raghu O and Philip J 2006 Thermal properties of paint coatings on different backings using a scanning photo acoustic technique Meas. Sci. Technol. 17 2945–9')
    lambda_paint = 1.38;
else
    lambda_paint = obj.HeatTransfer_params.HFS.lambda_paint;
end

if (strcmp(HFS.Type,'PCB'))
    if ~any(strcmp(hfs_data,'lambdax'))&&(obj.HeatTransfer_params.spatial_der)
        error('PIRT:Calculate_HeatTransfer: the lambdax must be introduced to compute the tangential heat transfer with the PCB model')
    end
    if ~any(strcmp(hfs_data,'lambday'))&&(obj.HeatTransfer_params.spatial_der)
        error('PIRT:Calculate_HeatTransfer: the lambday must be introduced to compute the tangential heat transfer with the PCB model')
    end
    PCB.lambdax = obj.HeatTransfer_params.HFS.lambdax;
    PCB.lambday = obj.HeatTransfer_params.HFS.lambday;
    k = [];
else
    if ~any(strcmp(hfs_data,'k'))
        error('PIRT:Calculate_HeatTransfer: the k must be introduced to compute the heat transfer')
    end
    k       = obj.HeatTransfer_params.HFS.k;
    PCB = [];
end

condition_data = fieldnames(Conditions);

if ~any(strcmp(condition_data,'V'))
    error('PIRT:Calculate_HeatTransfer: V must be introduced to compute the heat transfer')
end

if ~any(strcmp(condition_data,'I'))
    error('PIRT:Calculate_HeatTransfer: I must be introduced to compute the heat transfer')
end
if ~any(strcmp(condition_data,'L_char'))
    error('PIRT:Calculate_HeatTransfer: L must be introduced to compute the heat transfer')
end
if ~any(strcmp(condition_data,'Tamb'))
    error('PIRT:Calculate_HeatTransfer: Tamb must be introduced to compute the heat transfer')
end

if params.compute_St
    if ~any(strcmp(condition_data,'Uinf'))
        error('PIRT:Calculate_HeatTransfer: Uinf must be introduced to compute the heat transfer if the St number is desired')
    end
    if ~any(strcmp(condition_data,'rhoinf'))
        error('PIRT:Calculate_HeatTransfer: rhoinf must be introduced to compute the heat transfer if the St number is desired')
    end
end

s       = obj.HeatTransfer_params.HFS.s;
rho     = obj.HeatTransfer_params.HFS.rho;
cp      = obj.HeatTransfer_params.HFS.cp;
A       = obj.HeatTransfer_params.HFS.A;
epsilon = obj.HeatTransfer_params.HFS.epsilon;
L_char  = obj.HeatTransfer_params.conditions.L_char;
sigma   = obj.HeatTransfer_params.constants.sigma;
kp      = s_paint*lambda_paint;

params = fieldnames(obj.HeatTransfer_params);

if any(strcmp(params,'CustomQ'))
    customQ = obj.HeatTransfer_params.CustomQ;
else
    customQ = [];
end

end

function checkSizeMats(objective,Mat)

if length(size(objective)) == 3
    if length(size(Mat)) == 3
        if ~isequal(size(objective),size(Mat))
            error('PIRT:Calculate_HeatTransfer: The sizes of the images and the custom heat matrices must agree')
        end
    elseif length(size(Mat)) == 2
        if sum(size(Mat))==2
            if ~isfloat(Mat)
                error('PIRT:Calculate_HeatTransfer: Error in the customQ dimensions')
            end
        elseif ~isequal(size(objective,1:2),size(Mat))
            error('PIRT:Calculate_HeatTransfer: The sizes of the images and the custom heat matrices must agree')
        end
    else
        error('PIRT:Calculate_HeatTransfer: Error in the customQ dimensions')
    end
elseif length(size(objective)) == 2
    if length(size(Mat)) == 2
        if sum(size(Mat))==2
            if ~isfloat(Mat)
                error('PIRT:Calculate_HeatTransfer: Error in the customQ dimensions')
            end
        elseif ~isequal(size(objective,1:2),size(Mat))
            error('PIRT:Calculate_HeatTransfer: The sizes of the images and the custom heat matrices must agree')
        end
    else
        error('PIRT:Calculate_HeatTransfer: Error in the customQ dimensions')
    end
end


end

function [obj,q] = compute_tangential(q,Thot,obj,s,k,kp,PCB)


if isempty(obj.filter_params)
    % Compute derivatives with gradient
    fields = fieldnames(obj.HeatTransfer_params.conditions);
    if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
        % Obtain the derivative
        disp('----Computing spatial terms')
        [d2Tdx2_hot,d2Tdy2_hot,~] = Derivative_FD(Thot,1,0,obj.HeatTransfer_params.conditions.dx,obj.HeatTransfer_params.conditions.dy,0);
        disp('----Adding spatial terms')
        if isempty(PCB)
            q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
        else
            q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
        end

        if strcmp(obj.output.type,'file')
            disp('--> Saving d2Tdx2 and d2Tdy2 information into file')
            if endsWith(obj.output.path,'\')
                save(strcat([obj.output.path,'d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                save(strcat([obj.output.path,'d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
            else
                save(strcat([obj.output.path,'\d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                save(strcat([obj.output.path,'\d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
            end
            clear d2Tdx2_hot d2Tdy2_hot
        else
            obj.result.d2Tdx2_hot = d2Tdx2_hot;
            obj.result.d2Tdy2_hot = d2Tdy2_hot;
            clear d2Tdx2_hot d2Tdy2_hot
        end
    else
        warning('Tangencial-conduction could not be computed with gradient, the x and y resolution dx, dy are missing.')
    end
elseif ~isempty(find(obj.filter_params.filter==3))
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            if isfile(strcat([obj.output.path,'d2Tdx2.mat']))&&isfile(strcat([obj.output.path,'d2Tdy2.mat']))
                load(strcat([obj.output.path,'d2Tdx2.mat']))
                load(strcat([obj.output.path,'d2Tdy2.mat']))
                if isempty(PCB)
                    q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
                else
                    q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
                end
            else
                % Compute derivatives with gradient
                fields = fieldnames(obj.HeatTransfer_params.conditions);
                if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
                    warning('Tangencial-conduction term computed with gradient. sgolay32 filter data file was not found')
                    % Obtain the derivative
                    disp('----Computing spatial terms')
                    [d2Tdx2_hot,d2Tdy2_hot,~] = Derivative_FD(Thot,1,0,obj.HeatTransfer_params.conditions.dx,obj.HeatTransfer_params.conditions.dy,0);
                    disp('----Adding spatial terms')
                    if isempty(PCB)
                        q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
                    else
                        q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
                    end

                    if strcmp(obj.output.type,'file')
                        disp('--> Saving d2Tdx2 and d2Tdy2 information into file')
                        if endsWith(obj.output.path,'\')
                            save(strcat([obj.output.path,'d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                            save(strcat([obj.output.path,'d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                        else
                            save(strcat([obj.output.path,'\d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                            save(strcat([obj.output.path,'\d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                        end
                        clear d2Tdx2_hot d2Tdy2_hot
                    else
                        obj.result.d2Tdx2_hot = d2Tdx2_hot;
                        obj.result.d2Tdy2_hot = d2Tdy2_hot;
                        clear d2Tdx2_hot d2Tdy2_hot
                    end
                else
                    warning('Tangencial-conduction could not be computed with gradient, the x and y resolution dx, dy are missing.')
                end
            end
        else
            if isfile(strcat([obj.output.path,'\d2Tdx2.mat']))&&isfile(strcat([obj.output.path,'\d2Tdy2.mat']))
                load(strcat([obj.output.path,'\d2Tdx2.mat']))
                load(strcat([obj.output.path,'\d2Tdy2.mat']))
                if isempty(PCB)
                    q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
                else
                    q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
                end
            else
                % Compute derivatives with gradient
                fields = fieldnames(obj.HeatTransfer_params.conditions);
                if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
                    warning('Tangencial-conduction term computed with gradient. sgolay32 filter data file was not found')
                    % Obtain the derivative
                    disp('----Computing spatial terms')
                    [d2Tdx2_hot,d2Tdy2_hot,~] = Derivative_FD(Thot,1,0,obj.HeatTransfer_params.conditions.dx,obj.HeatTransfer_params.conditions.dy,0);
                    disp('----Adding spatial terms')
                    if isempty(PCB)
                        q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
                    else
                        q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
                    end

                    if strcmp(obj.output.type,'file')
                        disp('--> Saving d2Tdx2 and d2Tdy2 information into file')
                        if endsWith(obj.output.path,'\')
                            save(strcat([obj.output.path,'d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                            save(strcat([obj.output.path,'d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                        else
                            save(strcat([obj.output.path,'\d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                            save(strcat([obj.output.path,'\d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                        end
                        clear d2Tdx2_hot d2Tdy2_hot
                    else
                        obj.result.d2Tdx2_hot = d2Tdx2_hot;
                        obj.result.d2Tdy2_hot = d2Tdy2_hot;
                        clear d2Tdx2_hot d2Tdy2_hot
                    end
                else
                    warning('Tangencial-conduction could not be computed with gradient, the x and y resolution dx, dy are missing.')
                end
            end
        end

    else
        fields = fieldnames(obj.result);
        if and(any(strcmp(fields,'d2Tdx2_hot')),any(strcmp(fields,'d2Tdy2_hot')))
            disp('----Adding spatial terms')
            if isempty(PCB)
                q = q - (k+kp)*(obj.result.d2Tdx2_hot + obj.result.d2Tdy2_hot); % [W/m^2]
            else
                q = q - ((s*PCB.lambdax + kp)*obj.result.d2Tdx2_hot + (s*PCB.lambday + kp)*obj.result.d2Tdy2_hot);
            end
        else
            % Compute derivatives with gradient
            fields = fieldnames(obj.HeatTransfer_params.conditions);
            if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
                % Obtain the derivative
                disp('----Computing spatial terms')
                [d2Tdx2_hot,d2Tdy2_hot,~] = Derivative_FD(Thot,1,0,obj.HeatTransfer_params.conditions.dx,obj.HeatTransfer_params.conditions.dy,0);
                disp('----Adding spatial terms')
                if isempty(PCB)
                    q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
                else
                    q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
                end

                if strcmp(obj.output.type,'file')
                    disp('--> Saving d2Tdx2 and d2Tdy2 information into file')
                    if endsWith(obj.output.path,'\')
                        save(strcat([obj.output.path,'d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                        save(strcat([obj.output.path,'d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                    else
                        save(strcat([obj.output.path,'\d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                        save(strcat([obj.output.path,'\d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
                    end
                    clear d2Tdx2_hot d2Tdy2_hot
                else
                    obj.result.d2Tdx2_hot = d2Tdx2_hot;
                    obj.result.d2Tdy2_hot = d2Tdy2_hot;
                    clear d2Tdx2_hot d2Tdy2_hot
                end
            else
                warning('Tangencial-conduction could not be computed with gradient, the x and y resolution dx, dy are missing.')
            end
        end
    end
else
    % Compute derivatives with gradient
    fields = fieldnames(obj.HeatTransfer_params.conditions);
    if any(strcmp(fields,'dx'))&&any(strcmp(fields,'dy'))
        % Obtain the derivative
        disp('----Computing spatial terms')
        [d2Tdx2_hot,d2Tdy2_hot,~] = Derivative_FD(Thot,1,0,obj.HeatTransfer_params.conditions.dx,obj.HeatTransfer_params.conditions.dy,0);
        disp('----Adding spatial terms')
        if isempty(PCB)
            q = q - (k+kp)*(d2Tdx2_hot + d2Tdy2_hot); % [W/m^2]
        else
            q = q - ((s*PCB.lambdax + kp)*d2Tdx2_hot + (s*PCB.lambday + kp)*d2Tdy2_hot);
        end

        if strcmp(obj.output.type,'file')
            disp('--> Saving d2Tdx2 and d2Tdy2 information into file')
            if endsWith(obj.output.path,'\')
                save(strcat([obj.output.path,'d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                save(strcat([obj.output.path,'d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
            else
                save(strcat([obj.output.path,'\d2Tdx2.mat']),"d2Tdx2_hot",'-v7.3')
                save(strcat([obj.output.path,'\d2Tdy2.mat']),"d2Tdy2_hot",'-v7.3')
            end
            clear d2Tdx2_hot d2Tdy2_hot
        else
            obj.result.d2Tdx2_hot = d2Tdx2_hot;
            obj.result.d2Tdy2_hot = d2Tdy2_hot;
            clear d2Tdx2_hot d2Tdy2_hot
        end
    else
        warning('Tangencial-conduction could not be computed with gradient, the x and y resolution dx, dy are missing.')
    end
end


end

function [obj,q] = compute_unsteady(q,Thot,obj,rho,rho_p,s,s_p,cp,cp_p)

if isempty(obj.filter_params)
    % Compute derivatives with gradient
    fields = fieldnames(obj.HeatTransfer_params.conditions);
    if any(strcmp(fields,'dt'))
        % Obtain the derivative
        disp('----Computing unsteady terms')
        [~,~,dTdt_hot] = Derivative_FD(Thot,0,1,0,0,obj.HeatTransfer_params.conditions.dt);
        disp('----Adding unsteady terms')
        q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]

        if strcmp(obj.output.type,'file')
            disp('--> Saving dTdt information into file')
            if endsWith(obj.output.path,'\')
                save(strcat([obj.output.path,'dTdt.mat']),"dTdt_hot",'-v7.3')
            else
                save(strcat([obj.output.path,'\dTdt.mat']),"dTdt_hot",'-v7.3')
            end
            clear dTdt_hot
        else
            obj.result.dTdt_hot = dTdt_hot;
            clear dTdt_hot
        end
    else
        warning('Unsteady term could not be computed with gradient, the temporal resolution dt is missing.')
    end
elseif ~isempty(find(obj.filter_params.filter==3))
    if strcmp(obj.output.type,'file')
        if endsWith(obj.output.path,'\')
            if isfile(strcat([obj.output.path,'dTdt.mat']))
                load(strcat([obj.output.path,'dTdt.mat']))
                q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]
            else
                % Compute derivatives with gradient
                fields = fieldnames(obj.HeatTransfer_params.conditions);
                if any(strcmp(fields,'dt'))
                    warning('Unsteady term computed with gradient. sgolay32 filter data file was not found')
                    % Obtain the derivative
                    disp('----Computing unsteady terms')
                    [~,~,dTdt_hot] = Derivative_FD(Thot,0,1,0,0,obj.HeatTransfer_params.conditions.dt);
                    disp('----Adding unsteady terms')
                    q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]

                    if strcmp(obj.output.type,'file')
                        disp('--> Saving dTdt information into file')
                        if endsWith(obj.output.path,'\')
                            save(strcat([obj.output.path,'dTdt.mat']),"dTdt_hot",'-v7.3')
                        else
                            save(strcat([obj.output.path,'\dTdt.mat']),"dTdt_hot",'-v7.3')
                        end
                        clear dTdt_hot
                    else
                        obj.result.dTdt_hot = dTdt_hot;
                        clear dTdt_hot
                    end
                else
                    warning('Unsteady term could not be computed with gradient, the temporal resolution dt is missing.')
                end
            end
        else
            if isfile(strcat([obj.output.path,'\dTdt.mat']))
                load(strcat([obj.output.path,'\dTdt.mat']))
                q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]
            else
                % Compute derivatives with gradient
                fields = fieldnames(obj.HeatTransfer_params.conditions);
                if any(strcmp(fields,'dt'))
                    warning('Unsteady term computed with gradient. sgolay32 filter data file was not found')
                    % Obtain the derivative
                    disp('----Computing unsteady terms')
                    [~,~,dTdt_hot] = Derivative_FD(Thot,0,1,0,0,obj.HeatTransfer_params.conditions.dt);
                    disp('----Adding unsteady terms')
                    q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]

                    if strcmp(obj.output.type,'file')
                        disp('--> Saving dTdt information into file')
                        if endsWith(obj.output.path,'\')
                            save(strcat([obj.output.path,'dTdt.mat']),"dTdt_hot",'-v7.3')
                        else
                            save(strcat([obj.output.path,'\dTdt.mat']),"dTdt_hot",'-v7.3')
                        end
                        clear dTdt_hot
                    else
                        obj.result.dTdt_hot = dTdt_hot;
                        clear dTdt_hot
                    end
                else
                    warning('Unsteady term could not be computed with gradient, the temporal resolution dt is missing.')
                end
            end
        end

    else
        fields = fieldnames(obj.result);
        if any(strcmp(fields,'dTdt_hot'))
            disp('----Adding unsteady terms')
            q      = q - (rho*s*cp + rho_p*cp_p*s_p)*obj.result.dTdt_hot; % [W/m^2]
        else
            % Compute derivatives with gradient
            fields = fieldnames(obj.HeatTransfer_params.conditions);
            if any(strcmp(fields,'dt'))
                warning('Unsteady term computed with gradient. sgolay32 filter was not performed')
                % Obtain the derivative
                disp('----Computing unsteady terms')
                [~,~,dTdt_hot] = Derivative_FD(Thot,0,1,0,0,obj.HeatTransfer_params.conditions.dt);
                disp('----Adding unsteady terms')
                q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]

                if strcmp(obj.output.type,'file')
                    disp('--> Saving dTdt information into file')
                    if endsWith(obj.output.path,'\')
                        save(strcat([obj.output.path,'dTdt.mat']),"dTdt_hot",'-v7.3')
                    else
                        save(strcat([obj.output.path,'\dTdt.mat']),"dTdt_hot",'-v7.3')
                    end
                    clear dTdt_hot
                else
                    obj.result.dTdt_hot = dTdt_hot;
                    clear dTdt_hot
                end
            else
                warning('Unsteady term could not be computed with gradient, the temporal resolution dt is missing.')
            end
        end
    end
else
    % Compute derivatives with gradient
    fields = fieldnames(obj.HeatTransfer_params.conditions);
    if any(strcmp(fields,'dt'))
        warning('Unsteady term computed with gradient. sgolay32 filter was not performed')
        % Obtain the derivative
        disp('----Computing unsteady terms')
        [~,~,dTdt_hot] = Derivative_FD(Thot,0,1,0,0,obj.HeatTransfer_params.conditions.dt);
        disp('----Adding unsteady terms')
        q      = q - (rho*s*cp + rho_p*cp_p*s_p)*dTdt_hot; % [W/m^2]

        if strcmp(obj.output.type,'file')
            disp('--> Saving dTdt information into file')
            if endsWith(obj.output.path,'\')
                save(strcat([obj.output.path,'dTdt.mat']),"dTdt_hot",'-v7.3')
            else
                save(strcat([obj.output.path,'\dTdt.mat']),"dTdt_hot",'-v7.3')
            end
            clear dTdt_hot
        else
            obj.result.dTdt_hot = dTdt_hot;
            clear dTdt_hot
        end
    else
        warning('Unsteady term could not be computed with gradient, the temporal resolution dt is missing.')
    end
end
end
