function Nmod = findNmod(S,varargin)
% FINDNMOD Compute truncation for SVD modes.
%   FINDNMOD(S) ...
%   FINDNMOD(S,'Criterion',CRITERION,'Threshold',THRESHOLD) ... available
%   CRITERION:
%     - 'Spectrum': based on the energy spectrum of the SVD modes. This
%        quantity may reach asymptotically 1, for sufficiently large number
%        of modes. A reasonable THRESHOLD for the number of modes to retain
%        most valuable information can be set at 0.999 (default). This
%        CRITERION is formally equivalent to the variance fraction.
%     - 'Elbow': based on the ratio between forward and backward discrete
%        derivative of the reconstruction error. This quantity may reach
%        asymptotically 1, for sufficiently large number of modes. A
%        reasonable THRESHOLD for the number of modes to retain can be set
%        at 0.999 (default). This CRITERION is formally equivalent to look
%        for an elbow in the reconstruction error [1].
%     - 'HardThreshold': determining optimal turncation of Hard Threshold
%        for Matrix Denoising by Singular Values Hard Thresholding when
%        noise level is unknown. See optimal_SVHT_coef [2]
%
%   See also SVD, POD_FILTER, optimal_SVHT_coef
%
%   References:
%     [1] M. Raiola, S. Discetti and A. Ianiro, "On PIV random error 
%         minimization with optimal POD-based low order reconstruction",
%         https://e-archivo.uc3m.es/rest/api/core/bitstreams/ba43efba-8364-416d-9ae6-d7acaf936da1/content
%     [2] D. L. Donoho and M. Gavish, "The Optimal Hard Threshold for
%         Singular Values is 4/sqrt(3)", http://arxiv.org/abs/1305.5870
%
%   Author(s): M.Raiola, I.Robledo, R.Castellanos, S.Discetti, A.Ianiro
%   Copyright 2023 Universidad Carlos III de Madrid

[S,N,flag,threshold,beta] = parseinput(S,varargin{:});

lambda = S.^2;

switch flag
    case 1
        f = cumsum(lambda)./sum(lambda);
    case 2
        skip = 5;
        f = lambda(2:end)./lambda(1:end-1);
        f = smooth(f(skip:end-skip));  %Filtering
    case 3
        threshold = -optimal_SVHT_coef(beta,0) * median(S);
        f = -S;
end
Nmod = find(f>threshold,1);

if isempty(Nmod)
    Nmod = N;
end

fprintf('The suggested number of modes is: %d \n', Nmod)

end

%- Aux functions:
function [S,N,flag,threshold,beta] = parseinput(S,varargin)
beta = [];
threshold = [];

if isdiag(S)
    S = diag(S);
elseif ismatrix(S)
    error('PIRT:findNmod: The input S must be a diagonal matrix containing the eigenvalues')
else
    warning('PIRT:findNmod: The provided S is assumed as the vector of eigenvalues')
    S = S(:);
end
N = length(S);

if ~isempty(varargin)
    if any(strcmp(varargin,'Criterion'))
        idx = find(strcmp(varargin,'Criterion'));
        criterion = varargin{idx+1};
        if ~or(ischar(criterion),isstring(criterion))
            error('PIRT:findNmod: The introduced criterion must be text (string or char)')
        else
            switch criterion
                case {'Spectrum','spectrum',}
                    flag = 1;
                case {'Elbow','elbow'}
                    flag = 2;
                case {'HardThreshold','hardthreshold','hardThreshold'}
                    flag = 3;
                    if any(strcmp(varargin,'beta'))
                        idx = find(strcmp(varargin,'beta'));
                        beta = varargin{idx+1};
                        if ~isnumeric(beta)
                            error('PIRT:findNmod: The introduced beta must be numeric')
                        end
                    else
                        error('PIRT:findNmod: The data aspect ratio "beta" is required to use HardThreshold criterion')
                    end
                otherwise
                    error('PIRT:findNmod: The selected criterion is not valid')
            end
        end
    else
        warning('PIRT:findNmod: The Elbow criterion was selected as default')
        flag = 2;
    end

    if any(strcmp(varargin,'Threshold'))
        idx = find(strcmp(varargin,'Threshold'));
        threshold = varargin{idx+1};
        if ~isnumeric(threshold)
            error('PIRT:findNmod: The introduced Threshold must be numeric')
        end
        if threshold > 1
            warning('PIRT:findNmod: The introduced theshold can not be bigger than 1. The threshold selected was 0.999')
            threshold = 0.999;
        end
    elseif flag~=3
        warning('PIRT:findNmod: The default threshold is 0.999')
        threshold = 0.999;
    end

else
    warning('PIRT:findNmod: The Elbow criterion was selected as default with 0.999 as default threshold')
    flag = 2;
    threshold = 0.999;
end

end