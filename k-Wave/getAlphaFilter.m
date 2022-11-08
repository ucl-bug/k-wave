function alpha_filter = getAlphaFilter(kgrid, medium, filter_cutoff, taper_ratio)
%GETALPHAFILTER Create filter for medium.alpha_filter.
%
% DESCRIPTION:
%     getAlphaFilter uses getWin to create a Tukey window via rotation to
%     pass to the medium.alpha_filter input field of the first order
%     simulation functions (kspaceFirstOrder1D, kspaceFirstOrder2D, and
%     kspaceFirstOrder3D). This parameter is used to regularise time
%     reversal image reconstruction when absorption compensation is
%     included.
%       
% USAGE:
%     alpha_filter = getAlphaFilter(kgrid, medium, filter_cutoff)
%     alpha_filter = getAlphaFilter(kgrid, medium, filter_cutoff, taper_ratio)
% 
% INPUTS:
%     kgrid           - k-Wave grid object returned by kWaveGrid
%     medium          - k-Wave medium structure
%     filter_cutoff   - alpha_filter cutoff frequency [Hz], where
%                       filter_cutoff = [f_all] in 1D 
%                       filter_cutoff = [f_all] or [f_x, f_y] in 2D 
%                       filter_cutoff = [f_all] or [f_x, f_y, f_z] in 3D
%
%                       Any of the filter_cutoff inputs may be set to
%                       'max' to set the cutoff frequency to the maximum
%                       frequency supported by the grid
%
% OPTIONAL INPUTS
%     taper_ratio     - taper ratio for Tukey Window (default = 0.5)
%     
% OUTPUTS:
%     alpha_filter    - filter
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 6th May 2010
%     last update     - 4th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby
%
% See also getWin

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

% update the command line status
disp('Creating absorption filter...');

% check to see if the taper ratio is given
if nargin == 3
    taper_ratio = 0.5;
elseif nargin ~= 4
    error('Incorrect number of inputs.');
end

% update command line status
disp(['  taper ratio: ' num2str(taper_ratio)]);

% extract the dimensions
dim = numDim(kgrid.k);

% extract the alpha_filter cutoff frequencies
if numel(filter_cutoff) == 1
    filter_cutoff_x = filter_cutoff;
    if dim > 2
        filter_cutoff_z = filter_cutoff;
    end
    if dim > 1
        filter_cutoff_y = filter_cutoff;        
    end
elseif numel(filter_cutoff) ~= dim
    error(['Input filter_cutoff must have 1 or ' num2str(dim) ' elements for a ' num2str(dim) 'D grid.']);
else
    if dim == 1
        filter_cutoff_x = filter_cutoff{1};
    end
    if dim == 2
        filter_cutoff_x = filter_cutoff{1};
        filter_cutoff_y = filter_cutoff{2};        
    end    
    if dim == 3
        filter_cutoff_x = filter_cutoff{1};
        filter_cutoff_y = filter_cutoff{2}; 
        filter_cutoff_z = filter_cutoff{3};          
    end
end

% extract the maximium sound speed
c = max(medium.sound_speed(:));

% calculate the alpha_filter size in the z direction for 3D data
if dim > 2
    if strcmp(filter_cutoff_z, 'max')
        
        % set the alpha_filter size to be the same as the grid size
        filter_size_z = kgrid.Nz;
        filter_cutoff_z = kgrid.kz_max * c / (2 * pi);
        
    else
        
        % convert the cutoff frequency to a wavenumber
        k_cutoff_z = 2 * pi * filter_cutoff_z ./ c;
        
        % set the alpha_filter size
        filter_size_z = round(kgrid.Nz * k_cutoff_z / kgrid.kz(end));
        
        % check the alpha_filter size
        if filter_size_z > kgrid.Nz
            
            % set the alpha_filter size to be the same as the grid size
            filter_size_z = kgrid.Nz;   
            filter_cutoff_z = kgrid.kz_max * c / (2 * pi);
            
        end   
        
    end
end

% calculate the alpha_filter size in the y direction for 2 and 3D data
if dim > 1
    if strcmp(filter_cutoff_y, 'max')
        
        % set the alpha_filter size to be the same as the grid size
        filter_size_y = kgrid.Ny;
        filter_cutoff_y = kgrid.ky_max * c / (2 * pi);
        
    else
        
        % convert the cutoff frequency to a wavenumber
        k_cutoff_y = 2 * pi * filter_cutoff_y ./ c;
        
        % set the alpha_filter size
        filter_size_y = round(kgrid.Ny * k_cutoff_y / kgrid.ky(end));
        
        % check the alpha_filter size
        if filter_size_y > kgrid.Ny
            
            % set the alpha_filter size to be the same as the grid size
            filter_size_y = kgrid.Ny;
            filter_cutoff_y = kgrid.ky_max * c / (2 * pi);
            
        end
        
    end    
end

% calculate the alpha_filter size in the x direction for 1, 2, and 3D data
if strcmp(filter_cutoff_x, 'max')
    
    % set the alpha_filter size to be the same as the grid size
    filter_size_x = kgrid.Nx;
    filter_cutoff_x = kgrid.kx_max * c / (2 * pi);
    
else
    
    % convert the cutoff frequency to a wavenumber
    k_cutoff_x = 2 * pi * filter_cutoff_x ./ c;
    
    % set the alpha_filter size
    filter_size_x = round(kgrid.Nx * k_cutoff_x / kgrid.kx(end));
    
    % check the alpha_filter size
    if filter_size_x > kgrid.Nx
        
        % set the alpha_filter size to be the same as the grid size
        filter_size_x = kgrid.Nx;
        filter_cutoff_x = kgrid.kx_max * c / (2 * pi);
        
    end    
    
end  
    
% create the alpha_filter using getWin
switch dim
    case 1
        
        % create the alpha_filter
        filter_sec = getWin(filter_size_x, 'Tukey', 'Param', taper_ratio, 'Rotation', true);

        % enlarge the alpha_filter to the size of the grid
        alpha_filter  = zeros(kgrid.Nx, 1);
        x_index = round((kgrid.Nx - filter_size_x) / 2) + 1;
        alpha_filter(x_index:x_index + filter_size_x - 1) = filter_sec;  
        
        % update the command line status
        disp(['  filter cutoff: ' scaleSI(filter_cutoff_x) 'Hz']);
        
    case 2
        
        % create the alpha_filter
        filter_sec = getWin([filter_size_x, filter_size_y], 'Tukey', 'Param', taper_ratio, 'Rotation', true);

        % enlarge the alpha_filter to the size of the grid
        alpha_filter  = zeros(kgrid.Nx, kgrid.Ny);
        x_index = round((kgrid.Nx - filter_size_x) / 2) + 1;
        y_index = round((kgrid.Ny - filter_size_y) / 2) + 1;
        alpha_filter(x_index:x_index + filter_size_x - 1, y_index:y_index + filter_size_y - 1) = filter_sec;        
        
        % update the command line status
        disp(['  filter cutoff: ' scaleSI(filter_cutoff_x) 'Hz by ' scaleSI(filter_cutoff_y) 'Hz']);
        
    case 3
        
        % create the alpha_filter
        filter_sec = getWin([filter_size_x, filter_size_y, filter_size_z], 'Tukey', 'Param', taper_ratio, 'Rotation', true);

        % enlarge the alpha_filter to the size of the grid
        alpha_filter  = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
        x_index = round((kgrid.Nx - filter_size_x) / 2) + 1;
        y_index = round((kgrid.Ny - filter_size_y) / 2) + 1;
        z_index = round((kgrid.Nz - filter_size_z) / 2) + 1;
        alpha_filter(x_index:x_index + filter_size_x - 1, y_index:y_index + filter_size_y - 1, z_index:z_index + filter_size_z - 1) = filter_sec;
        
        % update the command line status
        disp(['  filter cutoff: ' scaleSI(filter_cutoff_x) 'Hz by ' scaleSI(filter_cutoff_y) 'Hz by ' scaleSI(filter_cutoff_z) 'Hz']);
        
end