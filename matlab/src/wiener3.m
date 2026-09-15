function [f,noise] = wiener3(varargin)
%WIENER3 3-D adaptive noise-removal filtering.
%   WIENER3 lowpass filters an intensity 3D images that has been degraded 
%   by constant power additive noise. WIENER3 uses a pixel-wise adaptive 
%   Wiener method based on statistics estimated from a local neighborhood 
%   of each pixel.
%
%   J = WIENER3(X,[M N L],NOISE) filters the 3D data X using pixel-wise 
%   adaptive Wiener filtering, using neighborhoods of size MxNxL 
%   (kernelsize) to estimate the local mean and standard deviation of the 
%   subdomain. If you omit the [M N L] argument, M N L default to 3. The
%   additive noise (Gaussian white noise) power is assumed to be NOISE.
%
%   [J,NOISE] = WIENER3(X,[M N,L]) also estimates the additive noise power
%   before doing the filtering. WIENER3 returns this estimate as NOISE.
%
%   Class Support
%   -------------
%   The input image I can be uint8, uint16, int16, double, or single.  The
%   output image J has the same class as I.
%
%   See also WIENER2, CONVN.
%
% Reference:
%   [1] Lim, Jae S. Two-Dimensional Signal and Image Processing, Englewood 
%       Cliffs, NJ, Prentice Hall, 1990, p. 548.
%
%   Author(s): R. Castellanos
%   Copyright 2023 Universidad Carlos III de Madrid


[g, nhood, noise] = ParseInputs(varargin{:});

classin = underlyingType(g);
classChanged = false;
if ~isUnderlyingType(g, 'double')
    classChanged = true;
    g = im2double(g);
end

% org = g;

% Estimate the local mean of f.
filter3 =@(stencil,X) convn(X,stencil,'same');

localMean = filter3(ones(nhood),g) / prod(nhood);

% Estimate of the local variance of f.
localVar = filter3(ones(nhood), g.^2) / prod(nhood) - localMean.^2;

% Estimate the noise power if necessary.
if (isempty(noise))
    noise = mean(localVar,'all');
end

% Compute result
% f = localMean + (max(0, localVar - noise) ./ ...
%           max(localVar, noise)) .* (g - localMean);
%
% Computation is split up to minimize use of memory
% for temp arrays.
f = g - localMean;
g = localVar - noise;
g = max(g, 0);
localVar = max(localVar, noise);
f = f ./ localVar;
f = f .* g;
f = f + localMean;
% aaaaaa = f;
% 
% time_series = reshape(org(351,200,:),[size(org,3),1]);
% Fs = 253;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = size(time_series,1);             % Length of signal
% t = (0:L-1)*T;        % Time vector
% f = Fs/L*(0:(L/2));
% 
% han_l    = 512;
% noverlap = han_l/2;
% P1 = pwelch(time_series,hanning(han_l),noverlap,f,Fs);
% semilogy(f,P1,'b')
% hold on
% 
% clear time_series
% time_series = reshape(aaaaaa(351,200,:),[size(aaaaaa,3),1]);
% Fs = 253;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = size(time_series,1);             % Length of signal
% t = (0:L-1)*T;        % Time vector
% f = Fs/L*(0:(L/2));
% 
% han_l    = 512;
% noverlap = han_l/2;
% P1 = pwelch(time_series,hanning(han_l),noverlap,f,Fs);
% semilogy(f,P1,'r')


if classChanged
    f = images.internal.changeClass(classin, f);
end


%%%
%%% Subfunction ParseInputs
%%%
function [g, nhood, noise] = ParseInputs(varargin)

nhood = [3 3 3];
noise = [];

switch nargin
    case 0
        error('wiener3: Not enought input arguments');
    case 1 % wiener3(I)
        g = varargin{1};
    case 2
        g = varargin{1};
        switch numel(varargin{2})
            case 1 % wiener3(I,noise)
                noise = varargin{2};
            case 3 % wiener3(I,[m n l])
                nhood = varargin{2};
            otherwise
                error('wiener3: Not valid kernel size. Either a scalar for uniform kernel size or 3-element vector');
        end
    case 3 % wiener3(I,[m n l],noise)
        g = varargin{1};
        nhood = varargin{2};
        noise = varargin{3};
    otherwise
        error('wiener3: Too many input arguments');
end