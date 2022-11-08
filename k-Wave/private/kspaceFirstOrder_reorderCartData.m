% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to reorder
%     the sensor points if a binary sensor mask is used for Cartesian
%     sensor mask nearest neighbour interpolation.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 2nd September 2012
%     last update - 21st March 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2019 Bradley Treeby

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

% update command line status
disp('  reordering Cartesian measurement data...');

if flags.record_p
    sensor_data.p = reorderBinarySensorData(sensor_data.p, reorder_index);
end

if flags.record_p_max
    sensor_data.p_max = reorderBinarySensorData(sensor_data.p_max, reorder_index);
end

if flags.record_p_min
    sensor_data.p_min = reorderBinarySensorData(sensor_data.p_min, reorder_index);
end

if flags.record_p_rms
    sensor_data.p_rms = reorderBinarySensorData(sensor_data.p_rms, reorder_index);
end

if flags.record_u
    
    % x-dimension
    sensor_data.ux = reorderBinarySensorData(sensor_data.ux, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy = reorderBinarySensorData(sensor_data.uy, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz = reorderBinarySensorData(sensor_data.uz, reorder_index);
    end
    
end

if flags.record_u_non_staggered
    
    % x-dimension
    sensor_data.ux_non_staggered = reorderBinarySensorData(sensor_data.ux_non_staggered, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_non_staggered = reorderBinarySensorData(sensor_data.uy_non_staggered, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_non_staggered = reorderBinarySensorData(sensor_data.uz_non_staggered, reorder_index);
    end
    
end

if flags.record_u_split_field
    
    % x-dimension
    sensor_data.ux_split_p = reorderBinarySensorData(sensor_data.ux_split_p, reorder_index);
    sensor_data.ux_split_s = reorderBinarySensorData(sensor_data.ux_split_s, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_split_p = reorderBinarySensorData(sensor_data.uy_split_p, reorder_index);
        sensor_data.uy_split_s = reorderBinarySensorData(sensor_data.uy_split_s, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_split_p = reorderBinarySensorData(sensor_data.uz_split_p, reorder_index);
        sensor_data.uz_split_s = reorderBinarySensorData(sensor_data.uz_split_s, reorder_index);
    end
    
end

if flags.record_u_max
    
    % x-dimension
    sensor_data.ux_max = reorderBinarySensorData(sensor_data.ux_max, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_max = reorderBinarySensorData(sensor_data.uy_max, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_max = reorderBinarySensorData(sensor_data.uz_max, reorder_index);
    end
    
end

if flags.record_u_min
    
    % x-dimension
    sensor_data.ux_min = reorderBinarySensorData(sensor_data.ux_min, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_min = reorderBinarySensorData(sensor_data.uy_min, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_min = reorderBinarySensorData(sensor_data.uz_min, reorder_index);
    end
    
end

if flags.record_u_rms
    
    % x-dimension
    sensor_data.ux_rms = reorderBinarySensorData(sensor_data.ux_rms, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_rms = reorderBinarySensorData(sensor_data.uy_rms, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_rms = reorderBinarySensorData(sensor_data.uz_rms, reorder_index);
    end
    
end

% note: intensity values are calculated before data re-ordering
if flags.record_I
    
    % x-dimension
    sensor_data.Ix = reorderBinarySensorData(sensor_data.Ix, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.Iy = reorderBinarySensorData(sensor_data.Iy, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.Iz = reorderBinarySensorData(sensor_data.Iz, reorder_index);
    end
    
end

if flags.record_I_avg
    
    % x-dimension
    sensor_data.Ix_avg = reorderBinarySensorData(sensor_data.Ix_avg, reorder_index);
    
    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.Iy_avg = reorderBinarySensorData(sensor_data.Iy_avg, reorder_index);
    end
    
    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.Iz_avg = reorderBinarySensorData(sensor_data.Iz_avg, reorder_index);
    end
    
end