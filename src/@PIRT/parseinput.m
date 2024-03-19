function obj = parseinput(obj,varargin)
% PARSEINPUT Reads and rearrange input from the user for PIRT class.
%   obj = PARSEINPUT('Tcold',A) assign the value of input A to the variable
%   Tcold. The input A must be numeric 2D or 3D matrix. The dimensions 1
%   and 2 (A(:,:,1), for example) correspond to the infrared snapshots
%   taken by the camera, being height x width of the image. Dimension 3
%   represents the time/sequence of snapshots.
%
%   obj = PARSEINPUT('Thot',B) assign the value of input B to the variable
%   Tcold. The input A must be numeric 2D or 3D matrix. Read previous
%   description.
%
%   ** PIRT requires 'Thot' and/or 'Tcold' to operate.
%
%   Examples:
%      % All possible inputs are considered.
%      A = [100 150 100; 150 100 150; 100; 150; 100];
%      B = [200 150 200; 150 200 150; 200; 150; 200];
%      Filter =
%
%      y = filloutliers(x, 'linear');

% It can have up to 7 extra inputs:
%       -Thot
%       -Tcold
%       -Filter
%       -Struct with the inputs associated with filtering
%       -CalculateHeatTransfer
%       -Struct with the inputs associated with calculating the
%            heat transfer number


% Checking for Wall temperature
if any(strcmp(varargin,'Thot'))
    idx = find(strcmp(varargin,'Thot'));
    if length(varargin)==idx
        error('PIRT:PIRT: After the identifier of Thot a matrix has to be introduced')
    end

    obj.Thot = varargin{idx+1};

    if ~isnumeric(obj.Thot)
        error('PIRT:PIRT:Thot must be a matrix')
    else
        if anynan(obj.Thot)
            warning('PIRT:PIRT: Thot contains NaN elements')
        end
        if any(isinf(obj.Thot))
            warning('PIRT:PIRT: Thot contains Inf elements')
        end
    end
end

% Checking for Adiabatic wall temperature
if any(strcmp(varargin,'Tcold'))
    idx = find(strcmp(varargin,'Tcold'));
    if length(varargin)==idx
        error('PIRT:PIRT: After the identifier of Tcold a matrix has to be introduced')
    end
    obj.Tcold = varargin{idx+1};
    if ~isnumeric(obj.Tcold)
        error('PIRT:PIRT:Tcold must be a matrix')
    else
        if anynan(obj.Tcold)
            warning('PIRT:PIRT: Tcold contains NaN elements')
        end
        if any(isinf(obj.Tcold))
            warning('PIRT:PIRT: Tcold contains Inf elements')
        end
    end
end

% Check if no input was provided
if and(isempty(obj.Thot),isempty(obj.Tcold))
    error('PIRT:PIRT: Input temperature is required. Either Thot or Tcold must be introduced to operate')
end

% Check for filter:
if any(strcmp(varargin,'Filter'))
    idx = find(strcmp(varargin,'Filter'));
    if length(varargin)==idx
        error('PIRT:PIRT: After the identifier of Filter a cell with the data has to be introduced')
    end

    filter_data = varargin{idx+1};
    input_data = obj.format_filter_inputs(filter_data);
    obj.filter_params = obj.parse_PIRT_filter(input_data{:});
    obj.Calculate_filter = true;
else
    warning('PIRT:PIRT:Filter No filter was selected')
    obj.filter_params = [];
end

cropping_points = [];

if any(strcmp(varargin,'Crop'))
    idx = find(strcmp(varargin,'Crop'));
    cropping_points = varargin{idx+1};
    if ~isfloat(cropping_points)
        error('PIRT:parse_PIRT_filter:The cropping points have to be introduces as a 2D')
    else
        if numel(cropping_points)~=4
            error('PIRT:parse_PIRT_filter:The cropping points introduced are the top left vertex and the bottom right vertex in a 2D matrix, being the first row the x coordinates and the bottom row the y coordinates')
    
        end
    end
end

obj.cropping_points = cropping_points;

if any(strcmp(varargin,'CalculateHeatTransfer'))
    if or(anynan(obj.Thot),anynan(obj.Tcold))
        error('PIRT:PIRT:In order to compute the Heat transfer both Thot and Tcold must be introduced')
    end

    idx = find(strcmp(varargin,'CalculateHeatTransfer'));
    if length(varargin)==idx
        error('PIRT:PIRT: After the identifier of CalculateHeatTransfer a cell array with the required inputs has to be introduced')
    end

    if any(strcmp(varargin,'HFS'))
        idx = find(strcmp(varargin,'HFS'));
        hfs_params = varargin{idx+1};
    else
        hfs_params = [];
    end

    if any(strcmp(varargin,'Conditions'))
        idx = find(strcmp(varargin,'Conditions'));
        conditions_params = varargin{idx+1};
    else
        conditions_params = [];
    end

    if any(strcmp(varargin,'Error'))
        idx = find(strcmp(varargin,'Error'));
        error_params = varargin{idx+1};
    else
        error_params = [];
    end

    input_data = obj.format_HeatTransfer_inputs(hfs_params,conditions_params,error_params);

    if any(strcmp(varargin,'TimeDer'))
        input_data = [input_data {'TimeDer'}];
    end
    if any(strcmp(varargin,'SpatialDer'))
        input_data = [input_data {'SpatialDer'}];
    end

    obj.HeatTransfer_params = obj.parse_PIRT_HeatTransfer(input_data{:});

    if any(strcmp(varargin,'Nu'))
        obj.HeatTransfer_params.compute_Nu = true;
    else
        obj.HeatTransfer_params.compute_Nu = false;
    end

    if any(strcmp(varargin,'St'))
        obj.HeatTransfer_params.compute_St = true;
    else
        obj.HeatTransfer_params.compute_St = false;
    end

    if any(strcmp(varargin,'h'))
        obj.HeatTransfer_params.compute_h = true;
    else
        obj.HeatTransfer_params.compute_h = false;
    end

    if any(strcmp(varargin,'CustomQ'))
        idx = find(strcmp(varargin,'CustomQ'));
        obj.HeatTransfer_params.CustomQ = varargin{idx+1};
        if ~iscell(obj.HeatTransfer_params.CustomQ)
            error('PIRT:PIRT: Custom heat tearms must be introduced in a cell')
        end
    end

    if and(and(~obj.HeatTransfer_params.compute_Nu,~obj.HeatTransfer_params.compute_St),~obj.HeatTransfer_params.compute_h)
        error('PIRT:PIRT: If the heat transfer option is selected either Nu, St or h must be introduced as an extra input')
    end

    obj.CalculateHeatTransfer = true;

    if any(strcmp(varargin,'CalculateHeatTransferError'))
        obj.CalculateHeatTransferError = true;
        idx = find(strcmp(varargin,'CalculateHeatTransferError'));
        if length(varargin)>idx
            errormethod = varargin{idx+1};
            switch errormethod
                case {'Moffat','moffat'}
                    obj.CalculateHeatTransferErrorMethod = 'Moffat';
                case {'Montecarlo','Monte Carlo','montecarlo','monte carlo'}
                    obj.CalculateHeatTransferErrorMethod = 'Montecarlo';
                otherwise
                    obj.CalculateHeatTransferErrorMethod = 'Moffat';
                    warning('PIRT:PIRT: Moffat method was selected to estimate the heat transfer error computation')
            end
        else
            obj.CalculateHeatTransferErrorMethod = 'direct';
            warning('PIRT:PIRT: direct method was selected to estimate the heat transfer error computation')
        end

    end
else
    disp('Heat transfer will not be computed')
    obj.HeatTransfer_params = NaN;
end