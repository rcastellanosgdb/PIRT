function parameters =  parse_PIRT_filter(obj,varargin)
%PARSE_PIRT_FILTER is a method that parses the inputs and generates the 
%   required flags and fields required for the execution of filtering in PIRT.
%   Required information: parse_PIRT_filter('Type',type)
%       Type: is a string or a vector of strings specifying the filters 
%           that will be performed
%   Optional information: parse_PIRT_filter('Parameters',Parameters,'Crop',Crop)
%       Parameters: is a cell or a vector of cells containing the relevant
%       parameters for each selected filter. The parameters for each filter 
%       must be introduced in the same order as the filters. Although not
%       required, most filters do require inputs so errors will arise
%       depending on the missing information.
%   This function also includes support functions for parsing the
%   parameters of each filter.
%   Returns:
%       parameters: struct with the information introduced in the format
%       that PIRT works in

if ~iscell(varargin)
    error('PIRT:parse_PIRT_filter:The input information must be introduced as a cell')
end

%% Check the filters

if ~any(strcmp(varargin,'Type'))
    error('PIRT:parse_PIRT_filter:Filters must be selected')
else
    idx = find(strcmp(varargin,'Type'));
    filters = varargin{idx+1};
end

if any(strcmp(varargin,'Parameters'))
    idx = find(strcmp(varargin,'Parameters'));
    params = varargin{idx+1};
else
    params = [];
end

if ~or(ischar(filters),isstring(filters))
    n_filters = length(filters);
    if n_filters~=length(params)
        error('PIRT:parse_PIRT_filter:If more than one filter is selected, the parameters have to be introduced as a cell with the same dimension as the number of filters')
    end

    POD_params = cell(n_filters,1);
    sgolay32_params = cell(n_filters,1);
    mPOD_params = cell(n_filters,1);
    wiener3_params = cell(n_filters,1);
    gaussian_params = cell(n_filters,1);
    cutoff_params = cell(n_filters,1);
    filter_type = zeros(n_filters,1);

    for i=1:n_filters
        [pod,golay,mpod,wiener,gaussian,cutoff,filter_type(i)] = select_filter(filters{i},params{i});
        POD_params{i} = pod;
        sgolay32_params{i} = golay;
        mPOD_params{i} = mpod;
        wiener3_params{i} = wiener;
        gaussian_params{i} = gaussian;
        cutoff_params{i} = cutoff;
    end

else
    [pod,golay,mpod,wiener,gaussian,cutoff,filter_type] = select_filter(filters,params);
    POD_params = cell(1,1);
    sgolay32_params = cell(1,1);
    mPOD_params = cell(1,1);
    wiener3_params = cell(1,1);
    gaussian_params = cell(1,1);
    cutoff_params = cell(1,1);
    POD_params{1} = pod;
    sgolay32_params{1} = golay;
    mPOD_params{1} = mpod;
    wiener3_params{1} = wiener;
    gaussian_params{1} = gaussian;
    cutoff_params{1} = cutoff;
end

%% Saving the parser information

parameters.filter = filter_type;
parameters.POD_params = POD_params;
parameters.sgolay32_params = sgolay32_params;
parameters.mPOD_params = mPOD_params;
parameters.gaussian_params = gaussian_params;
parameters.wiener3_params = wiener3_params;
parameters.cutoff_params = cutoff_params;

end

function [POD_params,sgolay32_params,mPOD_params,wiener3_params,gaussian_params...
    cutoff_params,filter_type] = select_filter(filter,params)

POD_params = [];
sgolay32_params = [];
mPOD_params = [];
wiener3_params = [];
gaussian_params = [];
cutoff_params = [];

if strcmp(filter,'POD')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the POD filter is selected, extra parameters are required')
    end
    POD_params = check_POD(params);
    filter_type = 1;
elseif strcmp(filter,'mPOD')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the m-POD filter is selected, extra parameters are required')
    end
    mPOD_params = checkmPOD(params);
    filter_type = 2;
elseif strcmp(filter,'sgolay32')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the sgolay32 filter is selected, extra parameters are required')
    end
    sgolay32_params = checksGolay(params);
    filter_type = 3;
elseif strcmp(filter,'wiener3')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the wiener3 filter is selected, extra parameters are required')
    end
    wiener3_params = checksWiener(params);
    filter_type = 4;
elseif strcmp(filter,'gaussian')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the gaussian filter is selected, extra parameters are required')
    end
    gaussian_params = checksGaussian(params);
    filter_type = 5;
elseif strcmp(filter,'cutoff')
    if ~iscell(params)
        error('PIRT:parse_PIRT_filter:If the Cutoff_Filter filter is selected, extra parameters are required')
    end
    cutoff_params = checksCutoff(params);
    filter_type = 6;
else
    error('PIRT:parse_PIRT_filter:A valid filter type must be selected')
end
end

function POD_params = check_POD(params)

if any(strcmp(params,'Threshold'))
idx = find(strcmp(params,'Threshold'));
POD_params.Threshold = params{idx+1};
    if ~isfloat(POD_params.Threshold)
        error('PIRT:parse_PIRT_filter:The Threshold introduced must be a float value')
    end
end
if any(strcmp(params,'Criterion'))
idx = find(strcmp(params,'Criterion'));
POD_params.Criterion = params{idx+1};
    if ~or(ischar(POD_params.Criterion),isstring(POD_params.Criterion))
        error('PIRT:parse_PIRT_filter:The Criterion introduced must be a string')
    end
end
if any(strcmp(params,'Nmod'))
idx = find(strcmp(params,'Nmod'));
POD_params.Nmod = params{idx+1};
    if ~isfloat(POD_params.Nmod)
        error('PIRT:parse_PIRT_filter:The Nmod introduced must be a float value')
    end
end
if any(strcmp(params,'beta'))
idx = find(strcmp(params,'beta'));
POD_params.beta = params{idx+1};
    if ~isfloat(POD_params.beta)
        error('PIRT:parse_PIRT_filter:TThe data aspect ratio "beta" must be a float value')
    end
end

if ~exist('POD_params','var')
    POD_params=struct([]);
end


end

function sgolay32_params = checksGolay(params)

if any(strcmp(params,'Kernel_size'))
    idx = find(strcmp(params,'Kernel_size'));
    sgolay32_params.w = params{idx+1};
    if ~isfloat(sgolay32_params.w)
        error('PIRT:parse_PIRT_filter:The Kernel_size introduced to the sgolay32 filter must be a float value or an array of float values')
    end
end

if any(strcmp(params,'h'))
    idx = find(strcmp(params,'h'));
    sgolay32_params.h = params{idx+1};
    if ~isfloat(sgolay32_params.h)
        error('PIRT:parse_PIRT_filter:The h (scaling) introduced to the sgolay32 filter must be a float value or an array of float values')
    end
end

if ~exist('sgolay32_params','var')
    sgolay32_params=struct([]);
end

end

function mPOD_params = checkmPOD(params)

if any(strcmp(params,'Type'))
    idx = find(strcmp(params,'Type'));
    switch params{idx+1}
        case {'PeakRem','Peak Removal','Peak removal'}
            if any(strcmp(params,'fpeaks'))
                idx = find(strcmp(params,'fpeaks'));
                mPOD_params.fpeak = params{idx+1};
                if ~isfloat(mPOD_params.fpeak)
                    error('PIRT:parse_PIRT_filter:The frequencies that are removed for the mPOD with peak remouval must be float values')
                end
            else
                error('PIRT:parse_PIRT_filter:The frequencies removed for the mPOD with peak remouval must be a introduced with label fpeaks')
            end
            if any(strcmp(params,'w'))
                idx = find(strcmp(params,'w'));
                mPOD_params.w = params{idx+1};
                if ~isfloat(mPOD_params.w)
                    error('PIRT:parse_PIRT_filter:The window size for the mPOD with peak remouval must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The window size for the mPOD with peak remouval must be a introduced with label w')
            end
            mPOD_params.type = 1;
        case {'Freq Deco','FreqDec','Frec Decoupling','Freq decoupling'}
            % if any(strcmp(params,'Threshold'))
            %     idx = find(strcmp(params,'Threshold'));
            %     mPOD_params.Threshold = params{idx+1};
            %     if ~isfloat(mPOD_params.Threshold)
            %         error('PIRT:parse_PIRT_filter:The Threshold for the mPOD filter must be a float value')
            %     end
            % end
            if any(strcmp(params,'N_regions'))
                idx = find(strcmp(params,'N_regions'));
                mPOD_params.Nf = params{idx+1};
                if ~isfloat(mPOD_params.Nf)
                    error('PIRT:parse_PIRT_filter:The number of regions for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The number of regions for the mPOD with frequency decoupling must be a introduced with label N_regions')
            end
            if any(strcmp(params,'f_min'))
                idx = find(strcmp(params,'f_min'));
                mPOD_params.fmin = params{idx+1};
                if ~isfloat(mPOD_params.fmin)
                    error('PIRT:parse_PIRT_filter:The minimum frequency for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The minimum for the mPOD with frequency decoupling must be a introduced with label f_min')
            end
            if any(strcmp(params,'f_max'))
                idx = find(strcmp(params,'f_max'));
                mPOD_params.fmax = params{idx+1};
                if ~isfloat(mPOD_params.fmax)
                    error('PIRT:parse_PIRT_filter:The maximum frequency for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The maximum for the mPOD with frequency decoupling must be a introduced with label f_max')
            end

            mPOD_params.type = 2;
        case {'Both','both'}

            if any(strcmp(params,'fpeaks'))
                idx = find(strcmp(params,'fpeaks'));
                mPOD_params.fpeak = params{idx+1};
                if ~isfloat(mPOD_params.fpeak)
                    error('PIRT:parse_PIRT_filter:The frequencies that are removed for the mPOD with peak remouval and frequency decoupling must be float values')
                end
            else
                error('PIRT:parse_PIRT_filter:The frequencies removed for the mPOD with peak remouval and frequency decoupling must be a introduced with label fpeaks')
            end

            if any(strcmp(params,'w'))
                idx = find(strcmp(params,'w'));
                mPOD_params.w = params{idx+1};
                if ~isfloat(mPOD_params.w)
                    error('PIRT:parse_PIRT_filter:The window size for the mPOD with peak remouval and frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The window size for the mPOD with peak remouval and frequency decoupling must be a introduced with label w')
            end

            % if any(strcmp(params,'Threshold'))
            %     idx = find(strcmp(params,'Threshold'));
            %     mPOD_params.Threshold = params{idx+1};
            %     if ~isfloat(mPOD_params.Threshold)
            %         error('PIRT:parse_PIRT_filter:The Threshold for the mPOD filter must be a float value')
            %     end
            % end
            if any(strcmp(params,'N_regions'))
                idx = find(strcmp(params,'N_regions'));
                mPOD_params.Nf = params{idx+1};
                if ~isfloat(mPOD_params.Nf)
                    error('PIRT:parse_PIRT_filter:The number of regions for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The number of regions for the mPOD with frequency decoupling must be a introduced with label N_regions')
            end
            if any(strcmp(params,'f_min'))
                idx = find(strcmp(params,'f_min'));
                mPOD_params.fmin = params{idx+1};
                if ~isfloat(mPOD_params.fmin)
                    error('PIRT:parse_PIRT_filter:The minimum frequency for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The minimum for the mPOD with frequency decoupling must be a introduced with label f_min')
            end
            if any(strcmp(params,'f_max'))
                idx = find(strcmp(params,'f_max'));
                mPOD_params.fmax = params{idx+1};
                if ~isfloat(mPOD_params.fmax)
                    error('PIRT:parse_PIRT_filter:The maximum frequency for the mPOD with frequency decoupling must be a float value')
                end
            else
                error('PIRT:parse_PIRT_filter:The maximum for the mPOD with frequency decoupling must be a introduced with label f_max')
            end
            mPOD_params.type = 3;
        otherwise
            error('PIRT:parse_PIRT_filter:If the m-POD filter is selected, Frequency decoupling, peak remouval or both must be selected')
    end
else
    error('PIRT:parse_PIRT_filter:The second input for the parameters in the m-POD must be its type')
end

if any(strcmp(params,'f_acquisition'))
    idx = find(strcmp(params,'f_acquisition'));
    mPOD_params.f_acquisition = params{idx+1};
    if ~isfloat(mPOD_params.f_acquisition)
        error('PIRT:parse_PIRT_filter:The aquisition frequency must be a float')
    end
else
    error('PIRT:parse_PIRT_filter:The aquisition frequency is required and its label is f_acquisition')
end

if any(strcmp(params,'Threshold'))
    idx = find(strcmp(params,'Threshold'));
    mPOD_params.Threshold = params{idx+1};
    if ~isfloat(mPOD_params.Threshold)
        error('PIRT:parse_PIRT_filter:The Threshold for the mPOD filter must be a float value')
    end
end

if ~exist('mPOD_params','var')
    mPOD_params=struct([]);
end

end

function wiener3_params = checksWiener(params)

if any(strcmp(params,'noise'))
    idx = find(strcmp(params,'noise'));
    wiener3_params.noise = params{idx+1};
    if ~isfloat(wiener3_params.noise)
        error('PIRT:parse_PIRT_filter:The noise introduced to the weiner3 filter must be a float value')
    end
end

if any(strcmp(params,'kernel'))
    idx = find(strcmp(params,'kernel'));
    wiener3_params.kernel = params{idx+1};
    if and(~isfloat(wiener3_params.kernel),length(wiener3_params.kernel)~=3)
        error('PIRT:parse_PIRT_filter:The kernel introduced to the weiner3 filter must be an array of 3 float values')
    end
end

if ~exist('wiener3_params','var')
    wiener3_params=struct([]);
end

end

function gaussian_params = checksGaussian(params)

if any(strcmp(params,'Sigma'))
    idx = find(strcmp(params,'Sigma'));
    gaussian_params.Sigma = params{idx+1};
    if ~isfloat(gaussian_params.Sigma)
        error('PIRT:parse_PIRT_filter:The standard deviation kernel introduced to the gaussian filter must be a float or a vector containing float values')
    end
    if (length(gaussian_params.Sigma)~=1)&&(length(gaussian_params.Sigma)~=3)
        error('PIRT:parse_PIRT_filter:The standard deviation kernel introduced to the gaussian filter must be only include 1 or 3 values')
    end
end

if any(strcmp(params,'FilterSize'))
    idx = find(strcmp(params,'FilterSize'));
    gaussian_params.FilterSize = params{idx+1};
    if ~isfloat(gaussian_params.FilterSize)
        error('PIRT:parse_PIRT_filter:The filter size introduced to the gaussian filter must be a float or a vector containing float values')
    end
    if (length(gaussian_params.FilterSize)~=1)&&(length(gaussian_params.FilterSize)~=3)
        error('PIRT:parse_PIRT_filter:The filter size introduced to the gaussian filter must be only include 1 or 3 values')
    end
end

if any(strcmp(params,'Padding'))
    idx = find(strcmp(params,'Padding'));
    switch params{idx+1}
        case 'replicate'
            gaussian_params.Padding = 'replicate';
        case 'circular'
            gaussian_params.Padding = 'circular';
        case 'symmetric'
            gaussian_params.Padding = 'symmetric';
        otherwise
            if isfloat(params{idx+1})
                gaussian_params.Padding = params{idx+1};
            else
                error('PIRT:parse_PIRT_filter: The introduced selection for the Padding in the imagaussfilt3 was not valid')
            end
    end
end

if any(strcmp(params,'FilterDomain'))
    idx = find(strcmp(params,'FilterDomain'));
    switch params{idx+1}
        case 'spatial'
            gaussian_params.FilterDomain = 'spatial';
        case 'frequency'
            gaussian_params.FilterDomain = 'frequency';
        case 'auto'
            gaussian_params.FilterDomain = 'auto';
        otherwise

            error('PIRT:parse_PIRT_filter: The introduced selection for the FilterDomain in the imagaussfilt3 was not valid')
    end
end

end

function cutoff_params = checksCutoff(params)

cutoff_params.Spatial = {};
cutoff_params.Temporal = {};

if ~any(strcmp(params,'Spatial'))&&~any(strcmp(params,'Temporal'))
    error('PIRT:parse_PIRT_filter: Spatial or Temporal fields must be introduced in the Cutoff filter')
end

if any(strcmp(params,'Spatial'))
    idx = find(strcmp(params,'Spatial'));
    spatial = params{idx+1};
    names = fieldnames(spatial);
    if any(strcmp(names,'fl'))
        if any(strcmp(names,'fh'))
            % Band-pass
            cutoff_params.Spatial{1} = spatial.fl;
            cutoff_params.Spatial{2} = spatial.fh;
        else
            % Low-pass
            cutoff_params.Spatial{1} = spatial.fl;
        end
    else
        if any(strcmp(names,'fh'))
            % High-pass
            cutoff_params.Spatial{1} = spatial.fh;
            cutoff_params.Spatial{2} = 'high';
        else
            warning('PIRT:parse_PIRT_filter: A spatial cutoff frequency was not specified, the default kernel size will be utilized')
        end
    end
end

if any(strcmp(params,'Temporal'))
    idx = find(strcmp(params,'Temporal'));
    temporal = params{idx+1};
    names = fieldnames(temporal);
    if ~any(strcmp(names,'fs'))
        error('PIRT:parse_PIRT_filter: The sampling frequency must be introduced for the temporal Cutoff filtering')
    end
    cutoff_params.Temporal{1} = temporal.fs;

    if any(strcmp(names,'fl'))
        if any(strcmp(names,'fh'))
            % Band-pass
            cutoff_params.Temporal{2} = temporal.fl;
            cutoff_params.Temporal{3} = temporal.fh;
        else
            % Low-pass
            cutoff_params.Temporal{2} = temporal.fl;
        end
    else
        if any(strcmp(names,'fh'))
            % High-pass
            cutoff_params.Temporal{2} = spatial.fh;
            cutoff_params.Temporal{3} = 'high';
        else
            error('PIRT:parse_PIRT_filter: A temporal cutoff frequency was not specified, either fl, fh or both must be introduced')
        end
    end
end

end
