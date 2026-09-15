% Copyright(C) 2019  UC3M - E.Díaz, R.Castellanos, C.Sanmiguel Vila, S.Discetti, A. Ianiro
%
% %******************************************************%
% %                                                      %
% %                       UC3M                           %
% %         Universidad Carlos III de Madrid             %
% %           Experimental Aerodynamics Lab              %
% %                                                      %
% %******************************************************%
%
% Author: E.Díaz
%
function [Tmap_new,Nmod] = multiscale_POD_filter(Tmap,params)

% Calculated according to "Multi-scale proper orthogonal decomposition of
% complex fluid flows". [Mendez et al. 2019 JFM].
Fs = params.f_acquisition;

[ny,nx,nt]  = size(Tmap);   % Dimensions of the problem
dt          = 1/Fs;  % Time step
df          = Fs/nt; % Frequency step

%% Reshape data and compute fluctuations
D       = reshape(Tmap,[],nt,1); %% Equation 2.1
clear Tmap;
D_mean  = mean(D,2); % Average
D       = D - D_mean; % Fluctuations

%% Compute time correlation matrix
K       = D'*D; % Section 2.1

%-- Compute Fourier Matrix (Section 2.2):
w_F     = exp(2*pi*sqrt(-1)/nt);
Fmat    = conj( w_F.^( ((1:nt)-1) .* ((1:nt)-1)' ) ./ sqrt(nt) ); % Fourier Matrix (Eq 2.13)

%-- Correlation matrix in Fourier space (Section 3.1):
K       = Fmat*K*Fmat; % (khat) Equation 3.1
clear w_F Fmat

%% Filter Bank (Generation of filter matrix)
switch params.type
    
    case 1
        H       = ones(nt); % Initialization of H matrix (frequency filter)
        w       = params.w; % Selection of window size
        f_vec   = params.fpeak; % Filtered frequencies
        % Remove the selected frequencies:
        for j = 1:length(f_vec)
            idx1 = round((f_vec(j)-w)/(df)); idx2 = round((f_vec(j)+w)/(df));
            H([idx1:idx2 nt-idx2:nt-idx1], :) = 0;
            H(:, [idx1:idx2 nt-idx2:nt-idx1]) = 0;
        end
        H       = sparse(H);
        clear w f_vec idx1 idx2 j filter

    case 2
        H       = zeros(nt); % Initialization of H matrix
        Nf      = params.Nf; fmin = params.fmin; fmax = params.fmax;
        idx0    = round(fmin/df);
        delta   = round(((fmax-fmin)/Nf)/df);
        for i = 1:Nf
            frange = [idx0+delta*(i-1) idx0+delta*i];
            H([frange(1)+1:frange(2) nt-frange(2)+1:nt-frange(1)],[frange(1)+1:frange(2) nt-frange(2)+1:nt-frange(1)]) = 1;
        end
        H       = sparse(H);
        clear Nf idx0 fmin fmax delta frange i filter
       
    case 3
        %-Frequency Decoupling
        H       = zeros(nt); % Initialization of H matrix
        Nf      = params.Nf; fmin = params.fmin; fmax = params.fmax;
        idx0    = round(fmin/df);
        delta   = round(((fmax-fmin)/Nf)/df);
        for i = 1:Nf
            frange = [idx0+delta*(i-1) idx0+delta*i];
            H([frange(1)+1:frange(2) nt-frange(2)+1:nt-frange(1)],[frange(1)+1:frange(2) nt-frange(2)+1:nt-frange(1)]) = 1;
        end
        %-Peak Removal
        w       = params.w; % Selection of window size
        f_vec   = params.fpeak; % Filtered frequencies
        for j = 1:length(f_vec)
            idx1 = round((f_vec(j)-w)/(df)); idx2 = round((f_vec(j)+w)/(df));
            H([idx1:idx2 nt-idx2:nt-idx1], :) = 0;
            H(:, [idx1:idx2 nt-idx2:nt-idx1]) = 0;
        end
        H       = sparse(H);
        clear w f_vec idx1 idx2 j idx0 fmin fmax delta frange i filter
        
    otherwise
        H       = ones(nt);
        clear filter
end


%% Filter and recover
K = K.*H; % Apply filter to the time correlation matrix in Fourier domain (frequency domain)
clear H

%-- Compute Fourier Matrix (Section 2.2):
w_F     = exp(2*pi*sqrt(-1)/nt);
Fmat    = conj( w_F.^( ((1:nt)-1) .* ((1:nt)-1)' ) ./ sqrt(nt) ); % Fourier Matrix (Eq 2.13)
clear w_F
% Get back to time domain with the inverse transform. K must be real, so
% take real() to get rid off the numerical complex residue.
%   - legacy_transform = false (default since v1.1): exact inverse,
%     K = real(F'*Khat*F').
%   - legacy_transform = true (behaviour of the code used in Robledo et al.
%     2025): the FORWARD matrix is applied again on both sides,
%     K = real(F*Khat*F). Since F*F is the index-reversal permutation this
%     returns the filtered correlation matrix with its time index reversed
%     (identical eigenvalues, time-reversed temporal modes), which damps the
%     reconstructed fluctuations (about -23% in Nu' for the sweeping-jet
%     dataset of the paper, see validation/REPORT.md).
if isfield(params,'legacy_transform') && params.legacy_transform
    K = real(Fmat*K*Fmat);   % legacy (paper) behaviour
else
    K = real(Fmat'*K*Fmat'); % exact inverse
end
clear Fmat

%-- Decompose the filtered time correlation matrix
[U,S,~] = svd(K, 'econ');
% Nmod = findNmod(S,'Criterion','HardThreshold','beta',2000/(399*604));
if isfield(params,'Threshold')
    Nmod = findNmod(S,'Criterion','Elbow','Threshold',params.Threshold);
else
    Nmod = findNmod(S,'Criterion','Elbow');
end
clear K

%-- Get spatial modes:
PHI = D*U/S.^0.5;
clear D
S = sparse(S);

%-- Reconstruct the matrix D:
Dnew = PHI(:,1:Nmod)*S(1:Nmod,1:Nmod).^0.5*U(:,1:Nmod)';
clear PHI Snew U


%-- Reshape the data matrix into the original form:
Tmap_new = reshape(Dnew + D_mean, [ny,nx,nt]);

% if Path
%     Tmap_new_mean = mean(Tmap_new,3);
%     savefast(strcat(Path,'\mPODdata\',Tcase,'_mPOD'),'Tmap_new','Tmap_new_mean');
% end

end
