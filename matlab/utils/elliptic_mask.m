function H = elliptic_mask(n,m,fn,fm,lowpass)
% ELLIPTIC_MASK common construction of the low-/high-pass spectral masks.
%   fn: normalised cut-off along the columns (m), fm: along the rows (n).
%   The semi-axes in frequency bins are floor(f*size/2). The grid is centred
%   on the DC bin of fftshift (index floor(size/2)+1), so the mask satisfies
%   H(k) = H(-k) and the filtered field is real.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

if fn>1 || fm>1
    error('elliptic_mask: The normalised cut-off frequencies can not exceed 1 (Nyquist)')
end
if fn<0 || fm<0
    error('elliptic_mask: The cut-off frequencies must be non-negative')
end

% Convert into actual sizes (frequency bins)
fn = floor(fn*m*0.5);
fm = floor(fm*n*0.5);
if fn<1 || fm<1
    warning('elliptic_mask: the cut-off is below the frequency resolution of the image, the mask degenerates')
end

% Frequency grid centred on the DC bin of fftshift
[kx,ky] = meshgrid((0:m-1)-floor(m/2),(0:n-1)-floor(n/2));

rho2 = (kx/fn).^2 + (ky/fm).^2;   % Inf where a semi-axis is 0
rho2(isnan(rho2)) = 0;            % 0/0 at the centre

if lowpass
    H = single(rho2<=1);
else
    H = single(rho2>=1);
end

end
