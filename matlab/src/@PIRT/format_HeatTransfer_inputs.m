function input_data = format_HeatTransfer_inputs(obj,HFS,Conditions,Error)
%FORMAT_HEATTRANSFER_INPUTS is a static method of class PIRT that translates the
%introduced structs containing the HT information to cells that will
%later be parsed
%   Inputs:
%       obj: required since it is a method
%       HFS: struct whith the foil sensor information. 
%       Conditions: struct whith the test conditions information. 
%       Error: struct with the error in heat transfer information
%   Returns:
%       input_data: cell with the translated information.

input_data = cell(1,4);
input_data{1} = 'HFS';
input_data{3} = 'Conditions';
input_data{5} = 'Error';

if ~isempty(HFS)
    hfs_values = struct2cell(HFS);
    hfs_names = fieldnames(HFS);

    hfs = cell(1,length(hfs_values)+length(hfs_names));

    hfs(1:2:end) = hfs_names;
    hfs(2:2:end) = hfs_values;

    input_data{2} = hfs;
else
    input_data{2} = [];
end

if ~isempty(Conditions)

    cond_values = struct2cell(Conditions);
    cond_names = fieldnames(Conditions);

    cond = cell(1,length(cond_values)+length(cond_names));

    cond(1:2:end) = cond_names;
    cond(2:2:end) = cond_values;

    input_data{4} = cond;
else
    input_data{4}=[];
end

if ~isempty(Error)

    err_values = struct2cell(Error);
    err_names = fieldnames(Error);

    err = cell(1,length(err_values)+length(err_values));

    err(1:2:end) = err_names;
    err(2:2:end) = err_values;

    input_data{6} = err;
else
    input_data=input_data(1:4);
end
end