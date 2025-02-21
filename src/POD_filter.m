function Xf = POD_filter(X,varargin)
% POD_FILTER Modal filtering based on Proper Orthogonal Decomposition.
%   POD_FILTER(X) ...
%
%   See also SVD, FINDNMOD.
%
% Reference:
%   [1] Sirovich, Lawrence (1987-10-01). "Turbulence and the dynamics of
%       coherent structures. I. Coherent structures". Quarterly of Applied
%       Mathematics. 45 (3): 561–571. doi:10.1090/qam/910462.
%
%   Author(s): I. Robledo, M. Raiola, R. Castellanos
%   Copyright 2023 Universidad Carlos III de Madrid

[Nmod,trunc_options] = parseinput(size(X),varargin{:});

[H,W,T] = size(X); % Extract dimensions of the temperature map
M       = reshape(X,[],size(X,3),1); clear X; % Rearrange the 3D matrix into 2D.
M_mean  = mean(M,2); % Average

% Substract the mean
X = M-M_mean; clear("M");

% Use the correlation matrix to minimize computational cost
if size(X,2)>=size(X,1) % More snapshots than temporal resolution
    [U,LAM] = eig(X*X','vector');
    % Sort eigenvalues and eigenvectors
    [lambda,ilam] = sort(LAM,'descend');
    U = U(:, ilam); % These are the spatial modes
    S = sqrt(diag(lambda)); clear LAM lambda ilam;
    % Calculate time coefficients
    V = X'*(U/S); clear X;

else % More spatial resolution than snapshots
    % Solve eigenvalue problem
    [V,LAM] = eig(X'*X,'vector');
    % Sort eigenvalues and eigenvectors
    [lambda,ilam] = sort(LAM,'descend');
    V = V(:,ilam); % These are the spatial modes
    S = sqrt(diag(lambda)); clear LAM lambda ilam;
    % Calculate time coefficients
    U = X*(V/S); clear X;
end

%-- Find the number of optimum modes to reconstruct
% Done if the user did not select Nmod as an input
if isempty(Nmod)
    if isempty(trunc_options)
        Nmod = findNmod(S);
    else
        Nmod = findNmod(S,trunc_options{:});
    end
end

% Reconstruct:
Xf = U(:,1:Nmod)*S(1:Nmod,1:Nmod)*V(:,1:Nmod)'; % Overwrite for memory issues
Xf  = reshape(Xf+M_mean, [H,W,T]); %Add back the mean and reshape
end

%- Aux func:
function [Nmod,trunc_options] = parseinput(images_shape,varargin)
trunc_options = {};
Nmod = [];
if ~isempty(varargin)
    if any(strcmp(varargin,'Nmod'))
        idx = find(strcmp(varargin,'Nmod'));
        Nmod = varargin{idx+1};
        if ~isnumeric(Nmod)
            error('PIRT:findNmod: The introduced Nmod must be numeric')
        end
    else
        trunc_options = varargin;
        if any(strcmp(varargin,'Criterion'))
            idx = find(strcmp(varargin,'Criterion'));
            criterion = varargin{idx+1};
            switch criterion
                case {'HardThreshold','hardthreshold','hardThreshold'}
                    if any(strcmp(varargin,'beta'))
                        idx = find(strcmp(varargin,'beta'));
                        beta = varargin{idx+1};
                        if ~isnumeric(beta)
                            error('PIRT:findNmod: The introduced beta must be numeric')
                        end
                    else
                        trunc_options{length(trunc_options)+1} = 'beta';
                        trunc_options{length(trunc_options)+1} = images_shape(3)/(images_shape(1)*images_shape(2));
                        warning('PIRT:findNmod: The data aspect ratio "beta" was not introduced and is selected as default the aspect ratio')
                    end
                otherwise
                    if any(strcmp(varargin,'Threshold'))
                        idx = find(strcmp(varargin,'Threshold'));
                        threshold = varargin{idx+1};
                        if ~isnumeric(threshold)
                            error('PIRT:findNmod: The introduced Threshold must be numeric')
                        end
                    end
            end
        end
    end
end
end