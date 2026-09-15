function C = sgolay32_coef(varargin)
% SGOLAY32_COEF Computes convolutional coefficient for Savitzky-Golay
% Filtering for 3D matrix and order 2 polynomial.
%    C = SGOLAY32_COEF(W) computes the 3D-convolution coefficients using a
%    framelength of W and a uniform spacing =1 in the three dimensions. The
%    order of the polynomial is fixed =2, so that it becomes a quadratic
%    filtering.
%    C = SGOLAY32_COEF(W,H) if H is scalar, considers a uniform scaling in
%    the three dimensions equals to H. If H is a 3-element vector,
%    considers non-uniform scaling.
%
%   See also SGOLAY, SGOLAYFILT, SGOLAY32_FILTER
%
%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.
%
%   Author(s): R. Castellanos, S. Discetti.
%   Copyright 2023 Universidad Carlos III de Madrid

[w,dx,dy,dz] = parseinput(varargin{:});

polynomial =@(i,j,k) [1,...
    i*dx, j*dy, k*dz,...
    i*j*dx*dy, i*k*dx*dz, j*k*dy*dz,...
    i*i*dx^2, j*j*dy^2, k*k*dz^2];

wcent = floor(w/2); %Centered region
rangex = -wcent(1):wcent(1);
rangey = -wcent(2):wcent(2);
rangez = -(-wcent(3):wcent(3));


M = zeros(numel(rangex)*numel(rangey)*numel(rangez),numel(polynomial(1,1,1)));

cont=0;
for k = rangez
    for j = rangey
        for i = rangex
            cont=cont+1;
            M(cont,:) = polynomial(i,j,k);
        end
    end
end

C  = (M'*M)\(M');
end

%- Aux:
function [w,dx,dy,dz] = parseinput(varargin)
switch nargin
    case 0
        w = 3*ones(1,3);
        h = ones(1,3);
    case 1 % SGOLAY32_COEF(W)
        w = varargin{1};
        h = ones(1,3);
    case 2 % SGOLAY32_COEF(W,H)
        w = varargin{1};
        h = varargin{2};
        if numel(h)==1
            h = repmat(h,[1,3]);
        elseif numel(h)~=3
            error('sgolay32_coef: Not valis spacing H. It must be a scalar (constant spacing) or a 3-element vector (non-uniform spacing)')
        end
    otherwise
        error('sgolay32_coef: Too many input arguments');
end

if numel(w)==1
    w = repmat(w,[1,3]);
elseif numel(w)~=3
    error('sgolay32_coef: Not valid kernel size W. Either a scalar for uniform kernel size or 3-element vector');
end

if any(mod(w,2)==0)
    error('SGOLAY32_FILTER: Kernel size must be odd');
end

dx = h(1); dy = h(2); dz = h(3);
end