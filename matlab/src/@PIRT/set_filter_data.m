function obj = set_filter_data(obj,varargin)
%SET_FILTER_DATA allows to insert further filter data for the later
%   calculations.
%
%   obj = SET_FILTER_DATA(Filter) saves the introduced information as
%   filter data to the later filtering. it must be introduced as a struct.
%
%   obj = SET_FILTER_DATA(Filter,'index',idx ) saves the introduced
%   information as filter data from the introduced index forward. It serves
%   to not eliminate previous introduced data.
%   obj = SET_FILTER_DATA(Filter,'append' ) appends the introduced filter
%   data.

if numel(varargin)==1
    if ~isstruct(varargin{1})
        error('PIRT:set_filter_data: The filter information must be introduced as a struct')
    end
    input_data = obj.format_filter_inputs(varargin{1});
    obj.filter_params = obj.parse_PIRT_filter(input_data{:});

else
    if ~isstruct(varargin{1})
        error('PIRT:set_filter_data: The filter information must be introduced as a struct in the first index')
    end

    if any(strcmp(varargin,'index'))
        if isempty(obj.filter_params)
            warning('PIRT:set_filter_data: No previous filter data was introduced so the introduced index will be ignored')
            index = 1;
        else
            idx = find(strcmp(varargin,'index'));
            index = varargin{idx+1};
            if floor(index)~=index
                error('PIRT:parse_PIRT_filter:The index introduced must be an integer')
            end
        end
        if ~isempty(new_filter_params.cropping_points)
            obj.filter_params.cropping_points = new_filter_params.cropping_points;
        end
    elseif any(strcmp(varargin,'append'))
        index = length(obj.filter_params.filter)+1;
    else
        warning('PIRT:set_filter_data: The extra information introduced was disregarded due to wrong inputs')
        index = 1;
    end

    input_data = obj.format_filter_inputs(varargin{1});
    new_filter_params = obj.parse_PIRT_filter(input_data{:});

    if index==1
        obj.filter_params.filter = new_filter_params.filter;
        obj.filter_params.POD_params = new_filter_params.POD_params;
        obj.filter_params.sgolay32_params = new_filter_params.sgolay32_params;
        obj.filter_params.mPOD_params = new_filter_params.mPOD_params;
        obj.filter_params.wiener3_params = new_filter_params.wiener3_params;
    else
        obj.filter_params.filter(index:index+length(new_filter_params.filter)-1) = new_filter_params.filter;
        obj.filter_params.POD_params(index:index+length(new_filter_params.filter)-1) = new_filter_params.POD_params;
        obj.filter_params.sgolay32_params(index:index+length(new_filter_params.filter)-1) = new_filter_params.sgolay32_params;
        obj.filter_params.mPOD_params(index:index+length(new_filter_params.filter)-1) = new_filter_params.mPOD_params;
        obj.filter_params.wiener3_params(index:index+length(new_filter_params.filter)-1) = new_filter_params.wiener3_params;
    end

end
end


