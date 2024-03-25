function Xf = Cutoff_Filter(X,varargin)
% CUTOFF_FILTER Performs a filtering operation based on a 2D low-pass,
% high-pass or a combination of both (band-pass) filters. The x direction
% in the filters is the row direction and the y the y direction
%   Cutoff_Filter(X): performs a low-pass filter operation with cutoff
%   frequencies equal to half of the allowed frequency range.
%   Cutoff_Filter(X,fc): performs a low-pass filter operation with cutoff
%   frequency fc in both x and y directions.
%   Cutoff_Filter(X,[fx,fy]): performs a low-pass filter operation with 
%   frequency fx in the x direction and fy in the y direction.
%   cutoff Cutoff_Filter(X,fc,'high'): performs a high-pass filter operation 
%   with cutoff frequency fc in both x and y directions.
%   Cutoff_Filter(X,[fx,fy],'high'): performs a high-pass filter operation 
%   with cutoff frequency fx in the x direction and fy in the y direction.
%   Cutoff_Filter(X,fl,fh): performs a band-pass filter  with cutoff
%   frequencies constant in x and y directions for the range [fl,fh].
%   Cutoff_Filter(X,[flx,fly],[fhx,fhy]): performs a band-pass filter 
%   with cutoff frequencies [flx,fly] for the lower cutoff range and
%   [fhx,fhy] for the higher cutoff range.
%   A single cutoff frequency for each range is also possible as
%   Cutoff_Filter(X,fl,[fhx,fhy]) or any other combination.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

% Obtain dimensions
dims = size(X);
if length(dims)==2
    n = dims(1);
    m = dims(2);
elseif length(dims)==3
    n = dims(1);
    m = dims(2);
    t = dims(3);
else
    error('Cutoff_Filter: only 2D and 3D matrix filtering supported')
end

if length(varargin)==1
    % Low-pass filter
    if ~isnumeric(varargin{1})
        error('Cutoff_Filter: the low-pass cutoff frequency must be numeric')
    end
    if length(varargin{1})==1
        H = lowpass_kernel(n,m,varargin{1},varargin{1});
    elseif length(varargin{1})==2
        H = lowpass_kernel(n,m,varargin{1}(2),varargin{1}(1));
    else
        error('Cutoff_Filter: the low-pass cutoff frequency can only be a scalar or a 2-element vector')
    end
elseif length(varargin)==2
    if strcmp(varargin{2},'high')
        % High-pass filter
        if ~isnumeric(varargin{1})
            error('Cutoff_Filter: the high-pass cutoff frequency must be numeric')
        end
        if length(varargin{1})==1
            H = highpass_kernel(n,m,varargin{1},varargin{1});
        elseif length(varargin{1})==2
            H = highpass_kernel(n,m,varargin{1}(2),varargin{1}(1));
        else
            error('Cutoff_Filter: the high-pass cutoff frequency can only be a scalar or a 2-element vector')
        end
    else
        % Band-pass filter
        if ~isnumeric(varargin{1})||~isnumeric(varargin{2})
            error('Cutoff_Filter: the low-pass and high-pass cutoff frequencies must be numeric')
        end
        % Obtain low-pass kernel
        if length(varargin{1})==1
            Hl = lowpass_kernel(n,m,varargin{1},varargin{1});
        elseif length(varargin{1})==2
            Hl = lowpass_kernel(n,m,varargin{1}(2),varargin{1}(1));
        else
            error('Cutoff_Filter: the low-pass cutoff frequency can only be a scalar or a 2-element vector')
        end
        % Obtain high-pass kernel
        if length(varargin{2})==1
            Hh = highpass_kernel(n,m,varargin{2},varargin{2});
        elseif length(varargin{2})==2
            Hh = highpass_kernel(n,m,varargin{2}(2),varargin{2}(1));
        else
            error('Cutoff_Filter: the high-pass cutoff frequency can only be a scalar or a 2-element vector')
        end
        H = Hl.*Hh;
        clear Hl Hh
    end
elseif isempty(varargin)
    H = lowpass_kernel(n,m,floor(n/4),floor(m/4));
else
    error('Cutoff_Filter: Inputs are not correctly introduced.')
end

% Obtain frequency content of the matrix
Freq = fft2(X);
% Center the frequency range
Freq = fftshift(Freq);

% figure()
% subplot(1,2,1)
% imagesc(H)
% axis equal
% subplot(1,2,2)
% imagesc(mean(abs(Freq),3))
% set(gca,'ColorScale','log')
% colorbar
% axis equal

% Apply the kernel
K = Freq.*H;

% Obtain the filtered matrix
Xf = abs(ifft2(ifftshift(K)));

end

