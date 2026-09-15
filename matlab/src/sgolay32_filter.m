function [Tfilt,dTdt,d2Tdx2,d2Tdy2] = sgolay32_filter(varargin)
%SGOLAY32_FILTER Savitzky-Golay filter for 3D matrix and order 2 polynomial
%   Tfilt = SGOLAY32_FILTER(X) smooths the data X using a Savitzky-Golay
%   (polynomial) smoothing filter.  A quadratic polynomial is considered
%   with uniform spacing =1 and framelength =3 in the three dimensions.
%
%   Tfilt= SGOLAY32_FILTER(X,KERNEL_SIZE) considers a uniform KERNEL_SIZE that
%   must be odd and greater than 2 (order of polynomial).  The size of the
%   input X in any of the 3 directions must be >= KERNEL_SIZE.
%
%   Note that if the polynomial order ORDER equals KERNEL_SIZE-1, no smoothing
%   will occur.
%
%   Tfilt= SGOLAY32_FILTER(X,KERNEL_SIZE,H) if H is scalar, considers a uniform
%   scaling in the three dimensions equals to H. If H is a 3-element vector
%   considers non-uniform scaling.
%
%   [Tfilt,dTdt,d2Tdx2,d2Tdy2] = SGOLAY32_FILTER(...) computes he time
%   derivative and the spatial second derivatives.
%
%   See also SGOLAY, SGOLAYFILT, SGOLAY32_COEF

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Castellanos, S. Discetti.
%   Copyright 2023 Universidad Carlos III de Madrid

[X,KERNEL_SIZE,dx,dy,dt] = parseinput(varargin{:});

C = sgolay32_coef(KERNEL_SIZE);

convcoef =@(col,w) reshape(col,[w(1) w(2) w(3)]);

expand = 1; % Hard coded to maintain the size of the original images

% Filtered temperature:
C0     = convcoef(C(1,:),KERNEL_SIZE);
Tfilt  = convn(X,C0,'same');
Tfilt = remove_edges( Tfilt ,KERNEL_SIZE,1);  % Remove edges and crop initial and final snapshots
% Time derivative:
Ct     = convcoef(C(4,:),KERNEL_SIZE);
dTdt   = convn(X,Ct,'same')./dt;
dTdt   = remove_edges(dTdt   ,KERNEL_SIZE,1);  % Remove edges and crop initial and final snapshots
% Second x derivative:
Cxx    = convcoef(C(9,:),KERNEL_SIZE);
d2Tdx2 = 2*convn(X,Cxx,'same')./dx^2;
d2Tdx2 = remove_edges(d2Tdx2 ,KERNEL_SIZE,1); % Remove edges and crop initial and final snapshots
% Second y derivative:
Cyy    = convcoef(C(8,:),KERNEL_SIZE);
d2Tdy2 = 2*convn(X,Cyy,'same')./dy^2;
d2Tdy2 = remove_edges(d2Tdy2 ,KERNEL_SIZE,1); % Remove edges and crop initial and final snapshots

end

%- Aux func:
function [X,w,dx,dy,dt] = parseinput(varargin)

switch nargin
    case 0
        error('SGOLAY32_FILTER: Not enough input arguments');
    case 1 %SGOLAY32_FILTER(X)
        w = 3;
        h = ones(1,3);
    case 2 % SGOLAY32_COEF(X,w)
        w = varargin{2};
        h = ones(1,3);
    case 3 % SGOLAY32_COEF(X,w,h)
        w = varargin{2};
        h = varargin{3};
        if numel(h)==1
            h = repmat(h,[1,3]);
        elseif numel(h)~=3
            error('SGOLAY32_FILTER: Not valis spacing H. It must be a scalar (constant spacing) or a 3-element vector (non-uniform spacing)')
        end
    otherwise
        error('SGOLAY32_FILTER: Too many input arguments');
end

X = varargin{1};
dx = h(1); dy = h(2); dt = h(3);

if numel(w)==1
    w = repmat(w,[1,3]);
elseif numel(w)~=3
    error('SGOLAY32_FILTER: Not valid kernel size W. Either a scalar for uniform kernel size or 3-element vector');
end

if any(mod(w,2)==0)
    error('SGOLAY32_FILTER: Kernel size must be odd');
end

end

function Xf = remove_edges(X,w,flag_crop_t)
win = floor(w/2); %Centered region

Xf = X;

%                          /|───────────────────/│
%                        /  |                 /  |
%                   1  /    |   Top         /    |
%                    /      |             / 2    |
%                  /        |           /        |
%                 │─────────|──────────│         |
%                 │         |          │         |
%                 │         |          │  Right  |
%                 │         |____________________|
%                 │        / Front     │        /
%                 │   3  /             │      /
%                 │    /               │    / 4
%                 │  /                 │  /
%                 │/___________________│/
% Walls of domain:
Xf(1:win(1),:,:)           = repmat(Xf(win(1)+1,:,:),win(1),1,1);   % Top
Xf((end-win(1)+1):end,:,:) = repmat(Xf(end-win(1),:,:),win(1),1,1); % Bottom
Xf(:,1:win(2),:)           = repmat(Xf(:,win(2)+1,:),1,win(2),1);   % Left
Xf(:,(end-win(2)+1):end,:) = repmat(Xf(:,end-win(2),:),1,win(2),1); % Right
Xf(:,:,1:win(3))           = repmat(Xf(:,:,win(3)+1),1,1,win(3));   % Front
Xf(:,:,(end-win(3)+1):end) = repmat(Xf(:,:,end-win(3)),1,1,win(3)); % Back

% Edges
Xf(1:win(1),1:win(2),:) = ...
    repmat(Xf(win(1)+1,win(2)+1,:),win(1),win(2),1);     % 1
Xf(1:win(1),(end-win(2)+1):end,:) = ...
    repmat(Xf(win(1)+1,end-win(2),:),win(1),win(2),1);   % 2
Xf((end-win(1)+1):end,1:win(2),:) = ...
    repmat(Xf(end-win(1),win(2)+1,:),win(1),win(2),1);   % 3
Xf((end-win(1)+1):end,(end-win(2)+1):end,:) = ...
    repmat(Xf(end-win(1),end-win(2),:),win(1),win(2),1); % 4

if flag_crop_t
    Xf = Xf(:,:,(win(3)+1):(end-win(3)));
end
end