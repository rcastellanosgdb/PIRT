function [asd] = Separate_filter(root,test,Nsnap,Ncold,mPODpoly,PODpoly)

folder=root;

Xrad=3;

if mPODpoly == 1
    hot_name_POD=strcat(folder,'\',test,'\postprocess\','Thot_m','POD');
    load(hot_name_POD);
    Mat = Thot_mPOD;
    clear Thot_mPOD
end
if PODpoly == 1
    hot_name_POD=strcat(folder,test,'\postprocess\','Thot_POD');
    load(hot_name_POD);
end

Matnew=Mat(:,:,1:Nsnap);
clear Mat

[varout,varout1,varout2,varout3,varout4,varout5,varout6]= Filtrospaziotempo(Matnew,Xrad);

if mPODpoly == 1
    hot_name_filter=strcat(folder,'\',test,'\postprocess\','Thot_m','filter');
end
if PODpoly == 1
    hot_name_filter=strcat(folder,'\',test,'\postprocess\','Thot_','filter');
end

savefast(hot_name_filter,'varout','varout1','varout2','varout3','varout4','varout5','varout6');%,'-v7.3')


%%%%%%%% Tcold

% cold_name=strcat(folder,'\',test,'\postprocess\','Tcold');
% 
% % cold_name=strcat(folder,'\postprocess\','Tcold_',test);
% 
% load(cold_name, 'Timage_cold')
% 
% Timage_coldnew = Timage_cold(:,:,1:Ncold);
% clear Timage_cold
% 
% [varout]= Filtrospaziotempo(Timage_coldnew,Xrad);
% 
% % cold_name_filter=strcat(folder,'\postprocess\','Tcold_',test,'filter');
% cold_name_filter=strcat(folder,'\',test,'\postprocess\','Tcold_','filter');
% 
% save(cold_name_filter,'varout','-v7.3')

asd = strcat('\n\tDone filtering ',test);

end