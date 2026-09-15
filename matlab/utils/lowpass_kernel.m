function H = lowpass_kernel(n,m,fn,fm)
% LOWPASS_KERNEL returns the kernel of a 2D low-pass filter with dimensions
% nxm with, cut-off frequencies of fn, fm. The kernel has an ellipsoidal
% form with the semi-major axes being fn (column direction, normalised to
% the Nyquist frequency) and fm (row direction). It is assumed that the
% frequency content has been shifted (fftshift) so that the zero frequency
% is at index floor(size/2)+1 in each dimension. n is the row direction and
% m the column direction.
%
%   v1.1: the frequency grid is now centred on the DC bin of fftshift and
%   the mask is exactly Hermitian-symmetric (previous versions used
%   (1:m)-m/2, offset by one bin, and mirrored one half of the mask).
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

H = elliptic_mask(n,m,fn,fm,true);

end
