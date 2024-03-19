% Trial Rodrigo Castellanos and Carlos Sanmiguel to validate the new
% Polyfilter code which is expected to be in the order to 10 times faster
% than the previous one:
clear; close all; clc;
rad = 3; % Number of adjacent poitns for the interpolation
load('H:\Airbus\RODRIGO\IR_Codes_src\Baseline\A1\postprocess\Thot')

Timage_hot=Timage_hot(:,:,1:30);
disp('New Polyfilter')
tic
[New.varout,New.varout1,New.varout2,New.varout3,...
    New.varout4,New.varout5,New.varout6] = PolyfilterCR(Timage_hot,rad);
toc

disp('Old Filtrospaziotempo')
tic
[Old.varout,Old.varout1,Old.varout2,Old.varout3,...
    Old.varout4,Old.varout5,Old.varout6] = Filtrospaziotempo(Timage_hot,rad);
toc
error.varout  = (New.varout  - Old.varout );
error.varout1 = (New.varout1 - Old.varout1);
% error.varout2 = (New.varout2 - Old.varout2)./Old.varout2;
% error.varout3 = (New.varout3 - Old.varout3)./Old.varout3;
% error.varout4 = (New.varout4 - Old.varout4)./Old.varout4;
% error.varout5 = (New.varout5 - Old.varout5)./Old.varout5;
% error.varout6 = (New.varout6 - Old.varout6)./Old.varout6;

% save('ComparePoly','New','Old','error');