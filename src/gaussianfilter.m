% Copyright(C) 2019  UC3M
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
function [asd] = gaussianfilter(root,test,Nsnap)

folder=root;

hot_name_POD=strcat(folder,'\postprocess\','Thot_',test,'POD');


load(hot_name_POD,'Mat');

Matnew=Mat(:,:,1:Nsnap);
clear Mat

Matnew=imgaussfilt3(Matnew);


hot_name_filter=strcat(folder,'\postprocess\','Thot_',test,'filter');


for k = 1:size(Matnew,3)
            if k==1
                varout3(:,:,k) = Matnew(:,:,k+1) - Matnew(:,:,k);          
            elseif k==size(Matnew,3)
                varout3(:,:,k) = Matnew(:,:,k) - Matnew(:,:,k-1); 
            else
                varout3(:,:,k) = (Matnew(:,:,k+1) - Matnew(:,:,k-1))/2;            
            end
end

% varout3 = diff(Matnew,1,3);
varout = Matnew;

save(hot_name_filter,'varout','varout3','-v7.3');


%%%%%%%% Tcold

cold_name=strcat(folder,'\postprocess\','Tcold_',test);

load(cold_name, 'Timage_cold')

Timage_coldnew = Timage_cold(:,:,1:Nsnap);
clear Timage_cold

% [varout]= Filtrospaziotempo(Timage_coldnew,Xrad);
[varout]= imgaussfilt3(Timage_coldnew);


cold_name_filter=strcat(folder,'\postprocess\','Tcold_',test,'filter');

save(cold_name_filter,'varout','-v7.3')

asd = strcat('Done filtering ',test);

end

