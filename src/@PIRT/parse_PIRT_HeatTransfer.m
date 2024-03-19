function parameters = parse_PIRT_HeatTransfer(obj,varargin)
%PARSE_PIRT_HEATTRANSFER is a method that parses the inputs and generates the 
% required flags and fields required for the calculation of the HT number.
%   Required information: parse_PIRT_HeatTransfer('HFS',HFS,'Conditions',Conditions)
%       HFS: is cell containing the introduced information of the thin foil
%       sensor, if any required inputs are missing an error will be given
%       Conditions: is cell containing the introduced information of the 
%       ambient conditions of the test, if any required inputs are missing
%       an error will be given.
%   Optional information: parse_PIRT_HeatTransfer('TimeDer','SpatialDer')
%       Error: cell containing the error information for the heat transfer
%       TimeDer: activates the time derivative computation
%       SpatialDer: activates the spatial derivative computation
%   Returns:
%       parameters: struct with the information introduced in the format
%       that PIRT works in

time_der = 0;
spatial_der = 0;

HFS = [];

conditions = [];

Error = [];

constants.g           = 9.80665; %[m/s^2] gravity constant
constants.sigma       = 5.67e-8;  % [W/(m^2 K^4)] Stefan-Boltzmann constant

%% Checking for time and spatial derivatives
if any(strcmp(varargin,'TimeDer'))
    time_der=1;
end
if any(strcmp(varargin,'SpatialDer'))
    spatial_der=1;
end

%% Obtaining HFS data

if any(strcmp(varargin,'HFS'))
    idx = find(strcmp(varargin,'HFS'));
    hfs_data = varargin{idx+1};
    clear HFS
    if ~isempty(hfs_data)
        if any(strcmp(hfs_data,'Type'))
            idx = find(strcmp(hfs_data,'Type'));
            HFS.Type = hfs_data{idx+1};
            if (~strcmp(HFS.Type,'Foil'))&&(~strcmp(HFS.Type,'PCB'))
                error('PIRT:parse_PIRT_HeatTransfer: The HTF Type has to be either Foil or PCB')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer: The HFS will be treated as a uniform thin foil')
            HFS.Type = 'Foil';
        end

        if any(strcmp(hfs_data,'lambdax'))
            idx = find(strcmp(hfs_data,'lambdax'));
            HFS.lambdax = hfs_data{idx+1};
            if ~isfloat(HFS.lambdax)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductance in the x direction must be introduced as a float value')
            end
        else
            if strcmp(HFS.Type,'PCB')&&spatial_der
                warning('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductance in the x direction must be introduced to calculate the tangential heat transfer in PCB and identified by "lambdax"')
            end
        end

        if any(strcmp(hfs_data,'lambday'))
            idx = find(strcmp(hfs_data,'lambday'));
            HFS.lambday = hfs_data{idx+1};
            if ~isfloat(HFS.lambday)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductance in the y direction must be introduced as a float value')
            end
        else
            if strcmp(HFS.Type,'PCB')&&spatial_der
                warning('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductance in the y direction must be introduced to calculate the tangential heat transfer in PCB and identified by "lambday"')
            end
        end

        if any(strcmp(hfs_data,'s'))
            idx = find(strcmp(hfs_data,'s'));
            HFS.s = hfs_data{idx+1};
            if ~isfloat(HFS.s)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS thickness "s" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The HFS thickness must be introduced to calculate the heat transfer and identified by "s"')
        end

        if any(strcmp(hfs_data,'cp'))
            idx = find(strcmp(hfs_data,'cp'));
            HFS.cp = hfs_data{idx+1};
            if ~isfloat(HFS.cp)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS "cp" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The HFS cp must be introduced to calculate the heat transfer and identified by "cp"')
        end

        if any(strcmp(hfs_data,'rho'))
            idx = find(strcmp(hfs_data,'rho'));
            HFS.rho = hfs_data{idx+1};
            if ~isfloat(HFS.rho)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS density "rho" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The HFS density must be introduced to calculate the heat transfer and identified by "rho"')
        end

        if and(any(strcmp(hfs_data,'H')),any(strcmp(hfs_data,'W')))
            idx = find(strcmp(hfs_data,'H'));
            HFS.H = hfs_data{idx+1};
            idx = find(strcmp(hfs_data,'W'));
            HFS.W = hfs_data{idx+1};
            if ~isfloat(HFS.H)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS height "H" must be introduced as a float value')
            end
            if ~isfloat(HFS.W)
                error('PIRT:parse_PIRT_HeatTransfer:The HFS width "W" must be introduced as a float value')
            end
            HFS.A = HFS.H*HFS.W;
        else
            if any(strcmp(hfs_data,'Area'))
                idx = find(strcmp(hfs_data,'Area'));
                HFS.A = hfs_data{idx+1};
                if ~isfloat(HFS.A)
                    error('PIRT:parse_PIRT_HeatTransfer:The HFS Area must be introduced as a float value')
                end

            else
                warning('PIRT:parse_PIRT_HeatTransfer:Either the HFS height AND width "H" and "W" or the HFS Area should be introduced to calculate the heat transfer')
            end
        end

        if any(strcmp(hfs_data,'epsilon'))
            idx = find(strcmp(hfs_data,'epsilon'));
            HFS.epsilon = hfs_data{idx+1};
            if ~isfloat(HFS.epsilon )
                error('PIRT:parse_PIRT_HeatTransfer:The HFS emissivity "epsilon" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The HFS emissivity must be introduced  to calculate the heat transfer and identified by "epsilon"')
        end

        if any(strcmp(hfs_data,'k'))
            idx = find(strcmp(hfs_data,'k'));
            HFS.k = hfs_data{idx+1};
            if ~isfloat(HFS.k )
                error('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductivity "k" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The HFS thermal conductivity must be introduced  to calculate the heat transfer and identified by "k"')
        end

%         HFS.k = HFS.s*HFS.lambda;
        HFS.m_plate = HFS.rho*HFS.s*HFS.A;
    else
        warning('PIRT:parse_PIRT_HeatTransfer:The HFS data must be introduced to calculate the heat transfer')
        HFS = [];
    end
else
    if isempty(obj.HeatTransfer_params.HFS)
        warning('PIRT:parse_PIRT_HeatTransfer:The HFS data must be introduced to calculate the heat transfer')
        HFS = [];
    end
end



%% Read conditions

if any(strcmp(varargin,'Conditions'))
    idx = find(strcmp(varargin,'Conditions'));
    cond = varargin{idx+1};
    clear conditions
    if ~isempty(cond)
        if any(strcmp(cond,'L'))
            idx = find(strcmp(cond,'L'));
            conditions.L_char = cond{idx+1};
            if ~isfloat(conditions.L_char)
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions characteristic length "L" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The Conditions characteristic length must be intoduced to calculate the heat transfer and identidied as "L"')
        end

        if any(strcmp(cond,'V'))
            idx = find(strcmp(cond,'V'));
            conditions.V = cond{idx+1};
            if ~isfloat(conditions.V)
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions voltage "V" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The Conditions voltage must be intoduced to calculate the heat transfer and identidied as "V"')
        end

        if any(strcmp(cond,'I'))
            idx = find(strcmp(cond,'I'));
            conditions.I = cond{idx+1};
            if ~isfloat(conditions.I)
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions current "I" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The Conditions current must intoduced to calculate the heat transfer and be identidied as "I"')
        end

        if any(strcmp(cond,'Uinf'))
            idx = find(strcmp(cond,'Uinf'));
            conditions.Uinf = cond{idx+1};
            if ~isfloat(conditions.Uinf)
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions free stream velocity "Uinf" must be introduced as a float value')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The Conditions free stream velocity must be intoduced to calculate the heat transfer if the St number is desired, it should be identified as "Uinf"')
        end

        if any(strcmp(cond,'Tamb'))
            idx = find(strcmp(cond,'Tamb'));
            conditions.Tamb = cond{idx+1};
            if ~isfloat(conditions.Tamb)
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions ambient temperature "Tamb" must be introduced as a 1D matrix')
            end
            if length(conditions.Tamb)~=2
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions parse_PIRT_HeatTransferambient temperature "Tamb" matrix must have 2 values; the hot and cold temperatures')
            end
            if or(~isfloat(conditions.Tamb(1)),~isfloat(conditions.Tamb(2)))
                error('PIRT:parse_PIRT_HeatTransfer:The Conditions ambient temperature "Tamb" values must be introduced as float')
            end
        else
            warning('PIRT:parse_PIRT_HeatTransfer:The Conditions ambient temperature must be identidied intoduced to calculate the heat transfer and as "Tamb"')
        end
    else
        warning('PIRT:parse_PIRT_HeatTransfer:The Conditions must be introduced to calculate the heat transfer')
        conditions = [];
    end
else
    if isempty(obj.HeatTransfer_params.conditions)
        warning('PIRT:parse_PIRT_HeatTransfer:The Conditions must be introduced to calculate the heat transfer')
        conditions = [];
    end
end

%% Reading error

if any(strcmp(varargin,'Error'))
    idx = find(strcmp(varargin,'Error'));
    err = varargin{idx+1};
    clear Error
    if ~isempty(err)
        if any(strcmp(err,'errorT'))
            idx = find(strcmp(err,'errorT'));
            Error.errorT = err{idx+1};
            if ~isfloat(Error.errorT)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the Temperatures must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the Temperatures must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorTamb'))
            idx = find(strcmp(err,'errorTamb'));
            Error.errorTamb = err{idx+1};
            if ~isfloat(Error.errorTamb)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the Tamb values must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the Tamb values must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorV'))
            idx = find(strcmp(err,'errorV'));
            Error.errorV = err{idx+1};
            if ~isfloat(Error.errorV)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the V must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the V value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorI'))
            idx = find(strcmp(err,'errorI'));
            Error.errorI = err{idx+1};
            if ~isfloat(Error.errorI)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the I must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the I value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorEpsilon'))
            idx = find(strcmp(err,'errorEpsilon'));
            Error.errorEpsilon = err{idx+1};
            if ~isfloat(Error.errorEpsilon)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the epsilon must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the epsilon value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorrho'))
            idx = find(strcmp(err,'errorrho'));
            Error.errorrho = err{idx+1};
            if ~isfloat(Error.errorrho)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the rho must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the rho value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorcp'))
            idx = find(strcmp(err,'errorcp'));
            Error.errorcp = err{idx+1};
            if ~isfloat(Error.errorcp)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the cp must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the cp value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errors'))
            idx = find(strcmp(err,'errors'));
            Error.errors = err{idx+1};
            if ~isfloat(Error.errors)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the s must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the s value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorA'))
            idx = find(strcmp(err,'errorA'));
            Error.errorAboard = err{idx+1};
            if ~isfloat(Error.errorAboard)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the Area of the sensor must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the Area of the sensor value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorkplate'))
            idx = find(strcmp(err,'errorkplate'));
            Error.errorkplate = err{idx+1};
            if ~isfloat(Error.errorkplate)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the k of the plate must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the k of the plate value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorLchar'))
            idx = find(strcmp(err,'errorLchar'));
            Error.errorLchar = err{idx+1};
            if ~isfloat(Error.errorLchar)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the characteristic L must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the characteristic L of the plate value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errork'))
            idx = find(strcmp(err,'errork'));
            Error.errork = err{idx+1};
            if ~isfloat(Error.errork)
                error('PIRT:parse_PIRT_HeatTransfer:The error in the k of the air must be introduced as a float value')
            end
        else
            error('PIRT:parse_PIRT_HeatTransfer:The error in the k of the air of value must be introduced to calculate the heat transfer error estimation')
        end

        if any(strcmp(err,'errorUinf'))
            idx = find(strcmp(err,'errorUinf'));
            Error.errorUinf = err{idx+1};
            if ~isfloat(Error.errorUinf)
                error('PIRT:parse_PIRT_HeatTransfer: The error in the Uinf must be introduced as a float value')
            end
        end

        if any(strcmp(err,'samples'))
            idx = find(strcmp(err,'samples'));
            Error.samples = err{idx+1};
            if rem(Error.samples,1)~=0
                error('PIRT:parse_PIRT_HeatTransfer: The number of samples for a Montecarlo error estimation must be introduced as an integer value')
            end
        end

        if time_der
            if any(strcmp(err,'errordTdt'))
                idx = find(strcmp(err,'errordTdt'));
                Error.errordTdt = err{idx+1};
                if ~isfloat(Error.errordTdt)
                    error('PIRT:parse_PIRT_HeatTransfer:The error in the dTdt must be introduced as a float value')
                end
            else
                error('PIRT:parse_PIRT_HeatTransfer:The error in the dTdt value must be introduced to calculate the heat transfer error estimation with the time derivative option')
            end
        end

        if spatial_der
            if any(strcmp(err,'errord2Tdx2'))
                idx = find(strcmp(err,'errord2Tdx2'));
                Error.errord2Tdx2 = err{idx+1};
                if ~isfloat(Error.errord2Tdx2)
                    error('PIRT:parse_PIRT_HeatTransfer:The error in the d2Tdx2 must be introduced as a float value')
                end
            else
                error('PIRT:parse_PIRT_HeatTransfer:The error in the d2Tdx2 value must be introduced to calculate the heat transfer error estimation with the spatial derivative option')
            end

            if any(strcmp(err,'errord2Tdy2'))
                idx = find(strcmp(err,'errord2Tdy2'));
                Error.errord2Tdy2 = err{idx+1};
                if ~isfloat(Error.errord2Tdy2)
                    error('PIRT:parse_PIRT_HeatTransfer:The error in the d2Tdy2 must be introduced as a float value')
                end
            else
                error('PIRT:parse_PIRT_HeatTransfer:The error in the d2Tdy2 value must be introduced to calculate the heat transfer error estimation with the spatial derivative option')
            end
        end

    else
        error('PIRT:parse_PIRT_HeatTransfer:The Conditions must be introduced to calculate the heat transfer')
        Error = [];
    end
end

%% Seting the parameters

parameters.time_der = time_der;
parameters.spatial_der = spatial_der;

parameters.HFS = HFS;

parameters.conditions = conditions;

parameters.constants = constants;

parameters.Error = Error;
end

