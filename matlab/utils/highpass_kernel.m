function H = highpass_kernel(n,m,fn,fm)
% HIGHPASS_KERNEL returns the kernel of a 2D high-pass filter with dimensions
% nxm with, cut-off frequencies of fn, fm. The kernel is the complement of
% the ellipsoidal low-pass kernel (see LOWPASS_KERNEL); the boundary of the
% ellipse is kept in both kernels.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

H = elliptic_mask(n,m,fn,fm,false);

end
