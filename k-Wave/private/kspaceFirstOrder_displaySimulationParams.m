% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to display
%     time steps and maximum supported frequency.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th July 2014
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2018 Bradley Treeby

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

% display time step information
disp(['  dt: ' scaleSI(kgrid.dt) 's, t_end: ' scaleSI(kgrid.t_array(end)) 's, time steps: ' num2str(kgrid.Nt)]);

% if using the elastic code, get the minimum sound speeds (not including
% zero if set for the shear speed)
if ~flags.elastic_code
    c_min = min(medium.sound_speed(:));
else
    c_min_comp = min(medium.sound_speed_compression(:));
    c_min_shear = min(medium.sound_speed_shear(medium.sound_speed_shear ~= 0));
end

% get suitable scaling factor
grid_size_metric = [kgrid.x_size, kgrid.y_size, kgrid.z_size];
[x_sc, scale, prefix] = scaleSI( min(grid_size_metric(grid_size_metric ~= 0)) ); %#ok<*ASGLU>
clear grid_size_metric;

% display the grid size and maximum supported frequency
switch kgrid.dim
    case 1
        
        % display grid size
        disp(['  input grid size: ' num2str(kgrid.Nx) ' grid points (' scaleSI(kgrid.x_size) 'm)']);
        
        % display maximum supported frequency
        disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max .* c_min ./ (2*pi) ) 'Hz']);         
    
    case 2
                
        % display grid size
        disp(['  input grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' grid points (' num2str(kgrid.x_size*scale) ' by ' num2str(kgrid.y_size*scale) prefix 'm)']);
        
        if ~flags.elastic_code
        
            % display maximum supported frequency
            if kgrid.kx_max == kgrid.ky_max
                disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max .* c_min ./ (2*pi) ) 'Hz']);
            else
                disp(['  maximum supported frequency: ' scaleSI( kgrid.kx_max .* c_min ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min ./ (2*pi) ) 'Hz']);
            end
            
        else

            % display the maximum supported frequency
            if kgrid.kx_max == kgrid.ky_max 
                disp(['  maximum supported compressional frequency: ' scaleSI( kgrid.k_max .* c_min_comp ./ (2*pi) ) 'Hz']);
                if isempty(c_min_shear)
                    disp('  maximum supported shear frequency: 0Hz');
                else
                    disp(['  maximum supported shear frequency: ' scaleSI( kgrid.k_max .* c_min_shear ./ (2*pi) ) 'Hz']);
                end
            else
                disp(['  maximum supported compressional frequency: ' scaleSI( kgrid.kx_max .* c_min_comp ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min_comp ./ (2*pi) ) 'Hz']);
                if isempty(c_min_shear)
                    disp('  maximum supported shear frequency: 0Hz');
                else
                    disp(['  maximum supported compressional frequency: ' scaleSI( kgrid.kx_max .* c_min_shear ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min_shear ./ (2*pi) ) 'Hz']);
                end
            end
            
        end

    case 3
        
        % display grid size
        disp(['  input grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' by ' num2str(kgrid.Nz) ' grid points (' num2str(kgrid.x_size*scale) ' by ' num2str(kgrid.y_size*scale) ' by ' num2str(kgrid.z_size*scale) prefix 'm)']); 
        
        if ~flags.elastic_code
        
            % display maximum supported frequency
            if (kgrid.kx_max == kgrid.kz_max) && (kgrid.kx_max == kgrid.ky_max)
                disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max .* c_min ./ (2*pi) ) 'Hz']);
            else
                disp(['  maximum supported frequency: ' scaleSI( kgrid.kx_max .* c_min ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.kz_max .* c_min ./ (2*pi) ) 'Hz']);
            end
        
        else
            
            % display the maximum supported frequency
            if (kgrid.kx_max == kgrid.kz_max) && (kgrid.kx_max == kgrid.ky_max)
                disp(['  maximum supported compressional frequency: ' scaleSI( kgrid.k_max .* c_min_comp ./ (2*pi) ) 'Hz']);
                if isempty(c_min_shear)
                    disp('  maximum supported shear frequency: 0Hz');
                else
                    disp(['  maximum supported shear frequency: ' scaleSI( kgrid.k_max .* c_min_shear ./ (2*pi) ) 'Hz']);
                end
            else
                disp(['  maximum supported frequency: ' scaleSI( kgrid.kx_max .* c_min_comp ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min_comp ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.kz_max .* c_min_comp ./ (2*pi) ) 'Hz']);
                if isempty(c_min_shear)
                    disp('  maximum supported shear frequency: 0Hz');
                else
                    disp(['  maximum supported frequency: ' scaleSI( kgrid.kx_max * c_min_shear ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.ky_max .* c_min_shear ./ (2*pi) ) 'Hz by ' scaleSI( kgrid.kz_max .* c_min_shear ./ (2*pi) ) 'Hz']);
                end
            end     
            
        end
        
end

% cleanup unused variables
clear c_min c_min_comp c_min_shear grid_size_metric prime_facs