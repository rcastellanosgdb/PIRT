function input_data = format_filter_inputs(obj,filter_data)
%FORMAT_FILTER_INPUTS is a static method of class PIRT that translates the
%introduced structs containing the filter information to cells that will
%later be parsed
%   Inputs:
%       obj: required since it is a method
%       filter_data: struct whith the filter information. Im must only have
%       2 fields: Type and Parameters where Type is the type of filter
%       selected and Parameters is a nested struct with the necessary
%       information. If filter data is a struct with several rows, the
%       translation will be done for each row independently
%   Returns:
%       input_data: cell with the translated information.


% Checking if the struct has more than 1 row
dims = size(filter_data);

%Seting the names
names = fieldnames(filter_data);
if length(names)==1
    names{2}='Parameters';
    filter_data.Parameters = struct([]);
end
input_data = cell(1,4);
input_data(1) = names(1);
input_data(3) = names(2);

values = struct2cell(filter_data);

if dims(2)>1
    data = cell(1,dims(2)); %Cell where the types will be stored
    parameters_data = cell(1,dims(2)); %Cell where the parameters will be stored
    for i=1:dims(2)
        current_data = values(:,:,i);
        if ischar(current_data{1})
            params = current_data{2};
            idx = 4;
            data(i) = current_data(1);
            idx_data = 2;
        else
            params = current_data{1};
            idx = 2;
            data(i) = current_data(2);
            idx_data = 4;
        end
        if isstruct(params)
            param_names = fieldnames(params);
            param_values = struct2cell(params);
            output = cell(1,length(param_names)+length(param_values)); %Cell where the information will be stored
            output(1:2:end) = param_names;
            output(2:2:end) = param_values;
            parameters_data{i} = output;
        else
            parameters_data{i} = {};
        end
    end
    input_data{idx} = parameters_data;
    input_data{idx_data} = data;
else
    %Checking which input is the params struct
    if isstruct(values{1})
        params = values{1};
        idx = 2;
        input_data(4) = values(2);
    else
        input_data(2) = values(1);
        idx  = 4;
        params = values{2};
    end
    param_names = fieldnames(params);
    param_values = struct2cell(params);
    output = cell(1,length(param_names)+length(param_values)); %Cell where the information will be stored
    output(1:2:end) = param_names;
    output(2:2:end) = param_values;
    input_data{idx} = output;
end
end