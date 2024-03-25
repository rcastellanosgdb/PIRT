function H = highpass_kernel(n,m,fn,fm)
% HIGHPASS_KERNEL returns the kernel of a 2D high-pass filter with dimensions
% nxm with, cut-off frequencies of fn, fm. The kernel has an ellipsoidal
% form with the sami-mayor axis being fn, fm. It is assumed that the
% frequency content has been shifted and centered around the 0 frequency
% content.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

% Check that the frequencies introduced are valid
if fn>floor(n/2)
    error('lowpass_kernel: The frequency cutoff for the y range is not valid')
end
if fm>floor(m/2)
    error('lowpass_kernel: The frequency cutoff for the x range is not valid')
end

% Create a meshgrid with the frequencies
[xx,yy] = meshgrid((1:m)-m/2,(1:n)-n/2);

% Find points inside frequency domain
xidx = xx>=fn*sqrt(1-yy.^2./fm^2);
if mod(m,2)==0
    xidx(:,1:floor(m/2)) = fliplr(xidx(:,floor(m/2)+1:end));
else
    xidx(:,1:floor(m/2)+1) = fliplr(xidx(:,floor(m/2)+1:end));
end
yidx = yy>=fm*sqrt(1-xx.^2./fn^2);
if mod(n,2)==0
    yidx(1:floor(n/2),:) = flip(yidx(floor(n/2)+1:end,:));
else
    yidx(1:floor(n/2)+1,:) = flip(yidx(floor(n/2)+1:end,:));
end

% Fill the kernel
H = single(xidx.*yidx);

end