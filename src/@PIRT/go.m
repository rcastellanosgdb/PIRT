function obj = go(obj)
%GO is a method that executes PIRT according to the inputs introduced

% Prepare the data
if length(size(obj.Tcold))>=3
    obj.Tcold = mean(obj.Tcold,3);
end

% Crop the images
if ~isempty(obj.cropping_points)
    points = obj.cropping_points;
    obj.Thot = obj.Thot(points(2,1):points(2,2),points(1,1):points(1,2),:);
    if ~isempty(obj.Tcold) 
        obj.Tcold = obj.Tcold(points(2,1):points(2,2),points(1,1):points(1,2),:);
    else
        warning('PIRT:PIRT: Only the hot images were introduced and cropped')
    end
end

if obj.Calculate_filter %Execute the filtering operation
    disp('*************************************************************')
    disp('******************** Temperature Filter *********************')
    disp('*************************************************************')
    
    obj.result.Thot_new = obj.Thot;

    obj.result.Nmod_hot = [];
    for i=1:length(obj.filter_params.filter)
        obj = obj.FilterTemperature(i);
    end
    if strcmp(obj.output.type,'file')
        disp('--> Saving information into file')
        Thot = obj.result.Thot_new;
        if endsWith(obj.output.path,'\')
            save(strcat([obj.output.path,'Thot_filtered.mat']),"Thot",'-v7.3')
        else
            save(strcat([obj.output.path,'\Thot_filtered.mat']),"Thot",'-v7.3')
        end
        clear Thot
        obj.result.Thot_new = [];
    end
else
    disp('******************** No Filter Applied *********************')
end

if obj.CalculateHeatTransfer %Execute the Nu computation
    disp('*************************************************************');
    disp('**************** Heat Transfer Calculation ******************');
    disp('*************************************************************');
    obj = obj.Calculate_HeatTransfer();

    if obj.CalculateHeatTransferError
        disp('*************************************************************');
        disp('************* Heat Transfer Error Calculation ***************');
        disp('*************************************************************');
        obj = obj.Calculate_HeatTransfer_Error();
        disp('-->DONE')
    end

end
end