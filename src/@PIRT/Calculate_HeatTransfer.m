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

[Thot,Tcold,s,rho,cp,k,A,epsilon,L_char,sigma,PCB,customQ] = parseinputs(obj);

disp('-- Calculating heat transfer maps');
% Check temperatures (convert to Kelvin if the average temperature is lower
% than an expected threshold):
threshold   = 100; % Rough estimation for minimum temperature in Celsius degrees:
Thot        = checkCelsius(Thot,threshold);
Tcold       = checkCelsius(Tcold,threshold);
Tamb_cold   = checkCelsius(obj.HeatTransfer_params.conditions.Tamb(1),threshold);
Tamb_hot    = checkCelsius(obj.HeatTransfer_params.conditions.Tamb(2),threshold);

ratio   = Tamb_hot/Tamb_cold; %[-]

%-- Heat Flux due to Joule Effect: qj = VI/A
qj      = obj.HeatTransfer_params.conditions.V*obj.HeatTransfer_params.conditions.I/(A); %[W/m^2]

%-- Radiative heat flux: qr = σ·ε·(Th⁴-Tc⁴)
qrad    = sigma.*epsilon.*(Thot.^4-Tamb_hot.^4); % [W/m^2]

%-- Internal energy
if obj.HeatTransfer_params.time_der
    if isempty(obj.filter_params)
        warning('Unsteady term not computed. sgolay32 filter has to be performed on the hot images')
        unsteady = 0;
    elseif ~isempty(find(obj.filter_params.filter==3))
        fields = fieldnames(obj.result);
        if any(strcmp(fields,'dTdt_hot'))
            disp('----Adding unsteady terms')
            unsteady = rho*s*cp*obj.result.dTdt_hot; % [W/m^2]
        else
            warning('Unsteady term not computed. sgolay32 filter has to be performed on the hot images')
            unsteady = 0;
        end
    else
        warning('Unsteady term not computed. sgolay32 filter has to be performed on the hot images')
        unsteady = 0;
    end
else
    unsteady=0;
end

%-- Tangencial conduction heat flux
if obj.HeatTransfer_params.spatial_der
    if isempty(obj.filter_params)
        warning('Tangencial-conduction term not computed. sgolay32 filter has to be performed on the hot images.')
        qk = 0;
    elseif ~isempty(find(obj.filter_params.filter==3))
        fields = fieldnames(obj.result);
        if and(any(strcmp(fields,'d2Tdx2_hot')),any(strcmp(fields,'d2Tdy2_hot')))
            disp('----Adding spatial terms')
            if isempty(PCB)
                qk = k*(obj.result.d2Tdx2_hot + obj.result.d2Tdy2_hot); % [W/m^2]
            else
                qk = s*(PCB.lambdax*obj.result.d2Tdx2_hot + PCB.lambday*obj.result.d2Tdy2_hot);
            end
        else
            warning('Tangencial-conduction term not computed. sgolay32 filter has to be performed on the hot images.')
            qk = 0;
        end
    else
        warning('Tangencial-conduction term not computed. sgolay32 filter has to be performed on the hot images.')
        qk = 0;
    end
else
    qk=0;
end

%Calculate heat
q = -qk - qrad + qj - unsteady;

if ~isempty(customQ)
    for i=1:length(customQ)
        checkSizeMats(Thot,customQ{i});
        q=q+customQ{i};
    end
end

% Calculate h (convective heat transfer coefficient) h [W/Km2]
h      = q./ ( Thot - (Tcold*ratio) ); %[W/K m^2]

if obj.HeatTransfer_params.compute_h
    obj.result.h      = h;
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
            Tfilm(:,:,nrange)   = (Thot(:,:,nrange) +Tcold_ref)/2; %[K]
        end
    catch
        Tfilm = (Thot+Tcold)/2;
    end
    % Air thermal conductivity as a function of T
    kair   = 1.5207E-11*Tfilm.^3-4.8574E-08*Tfilm.^2+1.0184E-04*Tfilm-3.9333E-04; clear Tfilm
    obj.result.Nu      = h.*L_char./kair;
end

if obj.HeatTransfer_params.compute_St
    Uinf    = obj.HeatTransfer_params.conditions.Uinf;
    obj.result.St      = h./(rho*cp*Uinf);
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

function [Thot,Tcold,s,rho,cp,k,A,epsilon,L_char,sigma,PCB,customQ] = parseinputs(obj)

HFS         = obj.HeatTransfer_params.HFS;
Conditions  = obj.HeatTransfer_params.conditions;
params      = obj.HeatTransfer_params;

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
    error('PIRT:Calculate_HeatTransfer: s must be introduced to compute the heat transfer')
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
end

s       = obj.HeatTransfer_params.HFS.s;
rho     = obj.HeatTransfer_params.HFS.rho;
cp      = obj.HeatTransfer_params.HFS.cp;
A       = obj.HeatTransfer_params.HFS.A;
epsilon = obj.HeatTransfer_params.HFS.epsilon;
L_char  = obj.HeatTransfer_params.conditions.L_char;
sigma   = obj.HeatTransfer_params.constants.sigma;

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