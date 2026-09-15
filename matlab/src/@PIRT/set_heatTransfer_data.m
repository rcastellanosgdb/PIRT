function obj = set_heatTransfer_data(obj,varargin)
%SET_HEATTRANSFER_DATA introduces to the object the introduced information to
%   calculate the heat transfer.
%   The available options are:
%       obj = obj.set_heatTransfer_data('HFS',HFS) which is a struct with the
%       heat foil sensor data
%       obj = obj.set_heatTransfer_data('Conditions',Conditions) which is a
%       struct with the conditions data
%       obj = obj.set_heatTransfer_data('TimeDer') which activates the time
%       derivative option for the HT calculation
%       obj = obj.set_heatTransfer_data('SpatialDer') which activates the
%       spatial derivatives
%       derivative option for the HT calculation
%   Any combination of the above mentioned inputs will also be admitted

N = 4;
combinations = dec2bin(0:2^N-1) - '0';

if any(strcmp(varargin,'TimeDer'))
    time_der=1;
else
    time_der=0;
end
if any(strcmp(varargin,'SpatialDer'))
    spatial_der=1;
else
    spatial_der=0;
end
if any(strcmp(varargin,'HFS'))
    HFS = 1;
    idx = find(strcmp(varargin,'HFS'));
    hfs_data = varargin{idx+1};
    input_data = obj.format_HeatTransfer_inputs(hfs_data,[],[]);
    hfs_data = input_data{2};
else
    HFS=0;
end
if any(strcmp(varargin,'Conditions'))
    Conditions = 1;
    idx = find(strcmp(varargin,'Conditions'));
    cond = varargin{idx+1};
    input_data = obj.format_HeatTransfer_inputs([],cond,[]);
    cond = input_data{4};
else
    Conditions = 0;
end

selection = [time_der spatial_der HFS Conditions];

idx = find(ismember(combinations,selection,'rows')==1);

switch idx
    case 1
        warning('PIRT:set_heatTransfer_data: No relevant heat transfer information was introduced')
    case 2
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond);
        obj.HeatTransfer_params.conditions = params.conditions;
    case 3
        params = obj.parse_PIRT_HeatTransfer('HFS',hfs_data);
        obj.HeatTransfer_params.HFS = params.HFS;
    case 4
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'HFS',hfs_data);
        obj.HeatTransfer_params.conditions = params.conditions;
        obj.HeatTransfer_params.HFS = params.HFS;
    case 5
        params = obj.parse_PIRT_HeatTransfer('SpatialDer');
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
    case 6
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'SpatialDer');
        obj.HeatTransfer_params.conditions = params.conditions;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
    case 7
        params = obj.parse_PIRT_HeatTransfer('HFS',hfs_data,'SpatialDer');
        obj.HeatTransfer_params.HFS = params.HFS;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
    case 8
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'HFS',hfs_data,'SpatialDer');
        obj.HeatTransfer_params.conditions = params.conditions;
        obj.HeatTransfer_params.HFS = params.HFS;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
    case 9
        params = obj.parse_PIRT_HeatTransfer('TimeDer');
        obj.HeatTransfer_params.time_der = params.time_der;
    case 10
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'TimeDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.conditions = params.conditions;
    case 11
        params = obj.parse_PIRT_HeatTransfer('HFS',hfs_data,'TimeDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.HFS = params.HFS;
    case 12
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'HFS',hfs_data,'TimeDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.HFS = params.HFS;
        obj.HeatTransfer_params.conditions = params.conditions;
    case 13
        params = obj.parse_PIRT_HeatTransfer('TimeDer','SpatialDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
    case 14
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'TimeDer','SpatialDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
        obj.HeatTransfer_params.conditions = params.conditions;
    case 15
        params = obj.parse_PIRT_HeatTransfer('HFS',hfs_data,'TimeDer','SpatialDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
        obj.HeatTransfer_params.HFS = params.HFS;
    case 16
        params = obj.parse_PIRT_HeatTransfer('Conditions',cond,'HFS',hfs_data,'TimeDer','SpatialDer');
        obj.HeatTransfer_params.time_der = params.time_der;
        obj.HeatTransfer_params.spatial_der = params.spatial_der;
        obj.HeatTransfer_params.HFS = params.HFS;
        obj.HeatTransfer_params.conditions = params.conditions;
end

