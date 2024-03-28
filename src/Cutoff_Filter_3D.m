function Xf = Cutoff_Filter_3D(X,varargin)
% CUTOFF_FILTER_3D performs a low/high/band-pass filtering operation in a
% 3D matrix in both the spatial (first and second dimension) and temporal
% dimension (third dim). It uses Cutoff_Filter as an auxiliary function to
% perform the spatial filtering and Temporal_Cutoff_Filter to perform the
% temporal filtering operation. The frequencies introduced must be in
% Hertz.
%   Cutoff_Filter_3D(X,'Spatial',fc): performs a low-pass spatial filter
%   operation with cutoff frequency fc
%   Cutoff_Filter_3D(X,'Spatial',fc,'high'): performs a high-pass spatial
%   filter operation with cutoff frequency fc
%   Cutoff_Filter_3D(X,'Spatial',fl,fh): performs a band-pass spatial
%   filter operation with low-cutoff frequency fl and high cutoff frequency
%   fh.
%       fc, fl and fh may be introduced as a scalar in which case the
%       kernel will be a circle of radius fc centered at 0 frequency. Or
%       they can be 2-element vectors with the sizes of the x and y axes of
%       an ellipsoid.
%
%   Cutoff_Filter_3D(X,'Temporal',fs,fc): performs a low-pass temporal filter
%   operation with cutoff frequency fc sampled at fs
%   Cutoff_Filter_3D(X,'Temporal',fs,fc,'high'): performs a high-pass temporal
%   filter operation with cutoff frequency fc sampled at fs
%   Cutoff_Filter_3D(X,'Temporal',fs,fl,fh): performs a band-pass temporal
%   filter operation with low-cutoff frequency fl and high cutoff frequency
%   fh sampled at fs.
%
%   Author(s): I. Robledo
%   Copyright 2023 Universidad Carlos III de Madrid

if any(strcmp(varargin,'Spatial'))
    idx_spat = find(strcmp(varargin,'Spatial'));
else
    idx_spat = -1;
end
if any(strcmp(varargin,'Temporal'))
    idx_temp = find(strcmp(varargin,'Temporal'));
else
    idx_temp = -1;
end

if idx_spat>-1&&idx_temp>-1
    % Both temporal and spatial filters
    if idx_spat>idx_temp
        spatial_idx = idx_spat+1:length(varargin);
        temporal_idx = idx_temp+1:idx_spat-1;
    else
        spatial_idx = idx_spat+1:idx_temp-1;
        temporal_idx = idx_temp+1:length(varargin);
    end
    % Spatial filter
    if length(spatial_idx)==1 % Low-pass
        Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)});
    elseif length(spatial_idx)==2
        if strcmp(varargin{spatial_idx(2)},'high') % high-pass
            Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)},'high');
        else
            % Band-pass
            Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)},varargin{spatial_idx(2)});
        end
    else
        error('Cutoff_Filter_3D: Information for the spatial filter not correctly introduced')
    end
    % Temporal filter
    if length(temporal_idx)==2 % Low-pass
        Xf  = Temporal_Cutoff_Filter(Xf,varargin{temporal_idx(1)},varargin{temporal_idx(2)});
    elseif length(temporal_idx)==3
        if strcmp(varargin{temporal_idx(3)},'high') % high-pass
            Xf  = Temporal_Cutoff_Filter(Xf,varargin{temporal_idx(1)},varargin{temporal_idx(2)},'high');
        else
        % Band-pass
            Xf  = Temporal_Cutoff_Filter(Xf,varargin{temporal_idx(1)},varargin{temporal_idx(2)},varargin{temporal_idx(3)});
        end
    else
        error('Cutoff_Filter_3D: Information for the temporal filter not correctly introduced')
    end
else
    if idx_spat>-1
        % Only Spatial filtering
        spatial_idx = idx_spat+1:length(varargin);
        
        if length(spatial_idx)==1 % Low-pass
            Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)});
        elseif length(spatial_idx)==2
            if strcmp(varargin{spatial_idx(2)},'high') % high-pass
                Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)},'high');
            else
                % Band-pass
                Xf = Spatial_Cutoff_Filter(X,varargin{spatial_idx(1)},varargin{spatial_idx(2)});
            end
        else
            error('Cutoff_Filter_3D: Information for the spatial filter not correctly introduced')
        end
    elseif idx_temp>-1
        % Only Temporal filtering
        temporal_idx = idx_temp+1:length(varargin);

        if length(temporal_idx)==2 % Low-pass
            Xf  = Temporal_Cutoff_Filter(X,varargin{temporal_idx(1)},varargin{temporal_idx(2)});
        elseif length(temporal_idx)==3
            if strcmp(varargin{temporal_idx(3)},'high') % high-pass
                Xf  = Temporal_Cutoff_Filter(X,varargin{temporal_idx(1)},varargin{temporal_idx(2)},'high');
            else
                % Band-pass
                Xf  = Temporal_Cutoff_Filter(X,varargin{temporal_idx(1)},varargin{temporal_idx(2)},varargin{temporal_idx(3)});
            end
        else
            error('Cutoff_Filter_3D: Information for the temporal filter not correctly introduced')
        end
    else
        error('Cutoff_Filter_3D: neither Spatial nor Temporal filtering were selected')
    end
end

end