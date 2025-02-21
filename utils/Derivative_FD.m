function [dTdx2,dTdy2,dTdt] = Derivative_FD(T,spatial,temporal,dx,dy,dt)
% DERIVATIVE_FD computes the derivatives of the matrix X utilizing second
% order finite differences. The row direction is treated as x, the column
% direction is treated as y and the third direction is the temporal
% dimension.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

dTdx2 = []; dTdy2 = []; dTdt = [];

if length(size(T))~=3&&length(size(T))~=2
    warning('The size of the matrix is incompatible with this method, no derivatives computed')
    spatial=0;temporal=0;
end
if length(size(T))==2&&temporal==1
    warning('No temporal derivative computed, the third dimention is the temporal dimension')
end

if spatial == 1
    if ~isnumeric(dx)||length(dx)~=1
        error('The discretization dx must be an unique number')
    end
    if ~isnumeric(dy)||length(dy)~=1
        error('The discretization dx must be an unique number')
    end
    if size(T,1)<4
        error('The number of rows is too small for the discretization implemented')
    end
    if size(T,2)<4
        error('The number of columns is too small for the discretization implemented')
    end
    dTdy2 = zeros(size(T));
    % Derivative using centered formula
    dTdy2(:,2:end-1,:) = (T(:,3:end,:)-2*T(:,2:end-1,:)+T(:,1:end-2,:))/dy^2;
    % Left margin with second-order forward approximation
    dTdy2(:,1,:) = (T(:,3,:) - 2*T(:,2,:) + T(:,1,:)) / dy^2;
    % Right margin with second-order backward approximation
    dTdy2(:,end,:) = (T(:,end-2,:) - 2*T(:,end-1,:) + T(:,end,:)) / dy^2;

    dTdx2 = zeros(size(T));
    % Derivative using centered formula
    dTdx2(2:end-1,:,:) = (T(3:end,:,:)-2*T(2:end-1,:,:)+T(1:end-2,:,:))/dx^2;
    % Left margin with second-order forward approximation
    dTdx2(1,:,:) = (T(3,:,:) - 2*T(2,:,:) + T(1,:,:)) / dx^2;
    % Right margin with second-order backward approximation
    dTdx2(end,:,:) = (T(end-2,:,:) - 2*T(end-1,:,:) + T(end,:,:)) / dx^2;
end

if temporal == 1
    % Temporal
    if ~isnumeric(dt)||length(dt)~=1
        error('The discretization dt must be an unique number')
    end
    if size(T,1)<3
        error('The number of rows is too small for the discretization implemented')
    end
    dTdt = zeros(size(T));
    % Derivative using centered formula
    dTdt(:,:,2:end-1) = (T(:,:,3:end)-T(:,:,1:end-2))/(2*dt);
    % Left margin with forward approximation
    dTdt(:,:,1) = (-3*T(:,:,1)+4*T(:,:,2)-T(:,:,3))/(2*dt);
    % Right margin woth backward approximation
    dTdt(:,:,end) = (3*T(:,:,end)-4*T(:,:,end-1)+T(:,:,end-2))/(2*dt);
end

end

