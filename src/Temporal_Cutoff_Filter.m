function Xf = Temporal_Cutoff_Filter(X,fs,varargin)
% Temporal_Cutoff_Filter performs a low-pass, high-pass or band-pass filter
% of a 3D matrix X taking into account that the third axis is the temporal
% domain.
%   Temporal_Cutoff_Filter(X,fs,fc): performs a low-pass filter on the 3D
%   matrix X acquired at a sampling rate fs [Hz] with a cutoff frequency of fc
%   [Hz].
%   Temporal_Cutoff_Filter(X,fs,fc,'high'): performs a high-pass filter on the 3D
%   matrix X acquired at a sampling rate fs [Hz] with a cutoff frequency of fc
%   [Hz].
%   Temporal_Cutoff_Filter(X,fs,fl,fc): performs a band-pass filter on the 3D
%   matrix X acquired at a sampling rate fs [Hz] with a low cutoff frequency
%   of fl [Hz] and a high cutoff frequency of fh [Hz].
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

if length(size(X))~=3
    error('Temporal_Cutoff_Filter: only valid for 3D matrices')
end
if ~isnumeric(X)
    error('Temporal_Cutoff_Filter: the matrix must be numeric')
end
if ~isnumeric(fs)
    error('Temporal_Cutoff_Filter: the sampling frequency must be numeric')
end
if length(fs)~=1
    error('Temporal_Cutoff_Filter: the sampling frequency must be one value')
end

% Obtain the dimensions
[n,m,l] = size(X);

% Reshape the matrix to utilize Matlab's built-in functions
Xr = reshape(permute(X,[3,2,1]),l,[]);

if length(varargin)==1
    %Low-pass
    if ~isnumeric(varargin{1})
        error('Temporal_Cutoff_Filter: the cutoff frequency for the low-pass filter must be numeric')
    end
    Xfr = lowpass(Xr,varargin{1},fs); % Perform the filtering
elseif length(varargin)==2
    if strcmp(varargin{2},'high')
        % high-pass
        if ~isnumeric(varargin{1})
            error('Temporal_Cutoff_Filter: the cutoff frequency for the high-pass filter must be numeric')
        end
        Xfr = highpass(Xr,varargin{1},fs); % Perform the filtering
    else
        % band-pass
        if ~isnumeric(varargin{1})
            error('Temporal_Cutoff_Filter: the cutoff frequency for the low-pass filter must be numeric')
        end
        if ~isnumeric(varargin{2})
            error('Temporal_Cutoff_Filter: the cutoff frequency for the high-pass filter must be numeric')
        end
        Xfr = bandpass(Xr,[varargin{1},varargin{2}],fs); % Perform the filtering
    end
else
    error('Temporal_Cutoff_Filter: Invalid number of inputs')
end
% Restore original shape
Xf = permute(reshape(Xfr,l,m,[]),[3,2,1]);
end