% Copyright(C) 2019  UC3M - R.Castellanos, C.Sanmiguel Vila, S.Discetti, A. Ianiro
%
% %******************************************************%
% %                                                      %
% %                       UC3M                           %
% %         Universidad Carlos III de Madrid             %
% %           Experimental Aerodynamics Lab              %
% %                                                      %
% %******************************************************%
%
%
function obj = FilterTemperature(obj,idx)
%FILTERTEMPERATURE is a method of class PIRT that dpending on the filter
%selected calls the appropiate functions to compute the new matrices.
%   Inputs:
%       obj: class object containing the required information such as the
%       images to filter, the types of filters and filter parameters
%       idx: index that specifies which filter is to be computed
%   Returns:
%       obj: the object with the new filtered images in the obj.result
%       attribute

Thot = obj.result.Thot_new;

switch obj.filter_params.filter(idx)
    % 0 - No filter
    % 1 - POD filter
    % 2 - Multi-scale POD filter
    % 3 - Sgolay
    % 4 - Wiener
    case 0
        disp('** No filter applied **');

    case 1 % POD filter
        disp('-- POD filter Hot images')
        if length(size(Thot))<3
            error('PIRT:FilterTemperature: A third dimension is required to apply the POD filter')
        end
        obj = iteratePOD(obj,Thot,idx);
        disp('--> Thot DONE')
    case 2 % Multi-scale
        disp('-- Multi-scale POD filter Hot images');
        if length(size(Thot))<3
            error('PIRT:FilterTemperature: A third dimension is required to apply the mPOD filter')
        end
        [obj.result.Thot_new,obj.result.Nmod_hot(idx)] = multiscale_POD_filter(Thot,obj.filter_params.mPOD_params{idx});
        disp('--> Thot DONE')

    case 3 % sgolay32 filter
        fields = fieldnames(obj.filter_params.sgolay32_params{idx});
        disp('-- sgolay32 filter Hot images');
        if ismember('w',fields)
            if ismember('h',fields)
                [obj.result.Thot_new,obj.result.dTdt_hot,obj.result.d2Tdx2_hot,obj.result.d2Tdy2_hot] = ...
                    sgolay32_filter(Thot,obj.filter_params.sgolay32_params{idx}.w,obj.filter_params.sgolay32_params{idx}.h);
            else
                [obj.result.Thot_new,obj.result.dTdt_hot,obj.result.d2Tdx2_hot,obj.result.d2Tdy2_hot] = ...
                    sgolay32_filter(Thot,obj.filter_params.sgolay32_params{idx}.w);
            end
        else
            [obj.result.Thot_new,obj.result.dTdt_hot,obj.result.d2Tdx2_hot,obj.result.d2Tdy2_hot] = ...
                sgolay32_filter(Thot);
        end
        disp('--> Thot DONE')

    case 4 %wiener3 filter
        disp('-- Wiener3 filter Hot images')
        obj = iterateWiener(obj,Thot,idx);
        disp('--> Thot DONE')
    otherwise
        error('Please select one of the options available for filters')
end

end

function obj = iteratePOD(obj,Tmat,idx)
fields = fieldnames(obj.filter_params.POD_params{idx});
inputs = cell(1,2*length(fields));
inputs(1:2:end)=fields;
inputs(2:2:end)=struct2cell(obj.filter_params.POD_params{idx});

[obj.result.Thot_new] = POD_filter(Tmat,inputs{:});

end

function obj = iterateWiener(obj,Tmat,idx)
fields = fieldnames(obj.filter_params.wiener3_params{idx});

if ismember('kernel',fields)
    if ismember('noise',fields)
        obj.result.Thot_new = wiener3(Tmat,...
            obj.filter_params.wiener3_params{idx}.kernel,...
            obj.filter_params.wiener3_params{idx}.noise);
    else
        [obj.result.Thot_new,obj.result.noise_hot] = wiener3(Tmat,...
            obj.filter_params.wiener3_params{idx}.kernel);

    end
else
    if ismember('noise',fields)
        obj.result.Thot_new = wiener3(Tmat,...
            obj.filter_params.wiener3_params{idx}.noise);
    else
        [obj.result.Thot_new,obj.result.noise_hot] = wiener3(Tmat);
    end
end

end
