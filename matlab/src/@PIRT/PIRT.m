classdef PIRT
    % PIRT Master class to operate with PIRT Toolbox
    %
    %   See also the routines in @PIRT

    %   Author: I. Robledo, J. Alfaro, R. Castellanos
    %   Copyright 2023 Universidad Carlos III de Madrid.

    properties
        Thot (:,:,:) {mustBeNumeric} = [] % Snapshot matrix with wall temperature
        Tcold (:,:,:) {mustBeNumeric} = [] % Snapshot matrix with adiabatic wall temperature
        Calculate_filter {mustBeNumericOrLogical} = false % Flag to filter the temperature maps.
        filter_params % Parameters for temperature filtering
        CalculateHeatTransfer {mustBeNumericOrLogical}  = false % Flag to calculate heat transfer.
        CalculateHeatTransferError {mustBeNumericOrLogical}  = false % Flag to calculate heat transfer error.
        CalculateHeatTransferErrorMethod %Method for error estimation
        HeatTransfer_params % Parameters for Nusselt number calculation
        output % Utilized for controlling how data saving is handled
        result = [] % Struct with output values
        cropping_points = [] % 2D matrix with the cropping points
    end

    methods
        obj = parseinput(obj,varargin)
        obj = set_filter_data(obj,varargin)
        obj = set_heatTransfer_data(obj,varargin)
        input_data = format_filter_inputs(obj,filter_data)
        input_data = format_HeatTransfer_inputs(obj,HFS,Conditions,Error)
        obj = go(obj)
        obj = FilterTemperature(obj,idx)
        obj = Calculate_HeatTransfer(obj)
        parameters =  parse_PIRT_filter(obj,varargin);
        parameters = parse_PIRT_HeatTransfer(obj,varargin);
        obj = Calculate_HeatTransfer_Error(obj);

        function obj = PIRT(varargin)
            % PIRT() Creates the object for the class
            %   obj = PIRT() generates a PIRT object with all the
            %   properties being empty. The input runflag is automatically
            %   set to 'false' to avoid executing PIRT.

            if ~isempty(varargin)
                obj = obj.parseinput(varargin{:});
                disp('PIRT: Information saved, run obj = obj.go() to perform calculations')
            else
                warning('PIRT: No inputs provided')
            end
        end
    end

end