clc; close all; clear
addpath("src/")
addpath("utils/")

path = './SJ_processing/';

[Timage_hot,Timage_cold,HFS,Conditions,f_ac,dx,dy] = load_case(path);

%% POD+gauss+Sgolay
clear Filter

% Create Filters
Filter(1).Type='POD';
Filter(1).Parameters.Criterion = 'HardThreshold';
Filter(1).Parameters.beta = 2000/(399*604);
Filter(2).Type='gaussian';
Filter(2).Parameters.FilterSize = [9,3,1];
Filter(2).Parameters.Sigma = [3,3,0.1];
Filter(3).Type = 'sgolay32';
Filter(3).Parameters.Kernel_size = [5,5,3];
Filter(3).Parameters.h = [dx,dy,1/f_ac];

% Create Object
output = PIRT('Thot',Timage_hot, ...  %Introduce images
              'Tcold',Timage_cold,... %Introduce images
              'CalculateHeatTransfer','Nu','h', ... % Params to compute
              'TimeDer','SpatialDer',... % Terms to consider
              'HFS',HFS,'Conditions',Conditions, ... %TFS parameters
              'Filter',Filter); % Filter Info
% Begin Computation
output = output.go();

% Extract results
Nu = output.result.Nu;
Num = mean(Nu,3);
Nuf = mean(abs(Nu-Num),3);

% Plot average and fluctuating maps
figure()
subplot(1,2,1)
imagesc(Num)
axis equal
axis off
title('Mean Nu Map')
subplot(1,2,2)
imagesc(Nuf)
axis equal
axis off
title('Fulctuating Nu Map')

function [Timage_hot,Timage_cold,HFS,Conditions, f_ac,dx,dy] = load_case(path)
% Load the resolution an test conditions
load(strcat([path,'Resolution.mat']));
load(strcat([path,'TestConditions.mat']));
f_ac = 1/dt;

% Save the conditions data
Conditions.V = V;
Conditions.I = I;
Conditions.Uinf = Uinf;
Conditions.Tamb = Tamb+273.15;
Conditions.dt = dt;

% Load the images
load(strcat([path,'Tcold.mat']));
load(strcat([path,'Thot.mat']));

% Save the TFS data
file = fopen(strcat([path,'CONF_DATA.STR']));

fgetl(file);
fgetl(file);
fgetl(file);
fgetl(file);
fgetl(file);

Conditions.L = sscanf(fgetl(file),'%e');

fgetl(file);

HFS.s = sscanf(fgetl(file),'%e');
HFS.rho = sscanf(fgetl(file),'%e');
HFS.cp = sscanf(fgetl(file),'%e');
HFS.k = sscanf(fgetl(file),'%e')*HFS.s;
HFS.H = sscanf(fgetl(file),'%e');
HFS.W = sscanf(fgetl(file),'%e');
HFS.epsilon = sscanf(fgetl(file),'%e');
HFS.sides   = 2;
HFS.s_paint = 42e-6;
HFS.cp_paint = 3061.5;
HFS.rho_paint = 1261.175;
HFS.lambda_paint = 1.38;


% Convert to meters
dx=1/dx/1000;
dy=1/dy/1000;

fprintf('Data loaded for case %s\n',path)

end