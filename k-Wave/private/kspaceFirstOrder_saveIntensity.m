% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to
%     calculate the acoustic intensity from the time varying acoustic
%     pressure and particle velocity recorded during the simulation. The
%     particle velocity is first temporally shifted forwards by dt/2 using
%     a Fourier interpolant so both variables are on the regular
%     (non-staggered) grid.
%
%     This function is called before the binary sensor data is reordered
%     for cuboid corners, so it works for both types of sensor mask.
%
%     If using cuboid corners with kspaceFirstOrder3DC/G, the sensor data
%     may be given as a structure array, i.e.,
%     sensor_data(n).ux_non_staggered. In this case, the calculation of
%     intensity is applied to each cuboid sensor mask separately.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 4th September 2013
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2018 Bradley Treeby

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

%#ok<*SAGROW>

% loop through the number of sensor masks (this can be > 1 if using cuboid
% corners with kspaceFirstOrder3DC/G)
for sensor_mask_index = 1:length(sensor_data)

    % shift the recorded particle velocity to the regular
    % (non-staggered) temporal grid
    ux_sgt     = fourierShift(sensor_data(sensor_mask_index).ux_non_staggered, 1/2);
    if kgrid.dim > 1
        uy_sgt = fourierShift(sensor_data(sensor_mask_index).uy_non_staggered, 1/2);
    end
    if kgrid.dim > 2        
        uz_sgt = fourierShift(sensor_data(sensor_mask_index).uz_non_staggered, 1/2);
    end 

    % compute the time varying intensity
    sensor_data(sensor_mask_index).Ix     = sensor_data(sensor_mask_index).p .* ux_sgt;
    if kgrid.dim > 1
        sensor_data(sensor_mask_index).Iy = sensor_data(sensor_mask_index).p .* uy_sgt;
    end
    if kgrid.dim > 2
        sensor_data(sensor_mask_index).Iz = sensor_data(sensor_mask_index).p .* uz_sgt;
    end  

    % calculate the time average of the intensity if required using the last
    % dimension (this works for both linear and cuboid sensor masks)
    if flags.record_I_avg
        sensor_data(sensor_mask_index).Ix_avg     = mean(sensor_data(sensor_mask_index).Ix, ndims(sensor_data(sensor_mask_index).Ix));
        if kgrid.dim > 1
            sensor_data(sensor_mask_index).Iy_avg = mean(sensor_data(sensor_mask_index).Iy, ndims(sensor_data(sensor_mask_index).Iy));
        end
        if kgrid.dim > 2
            sensor_data(sensor_mask_index).Iz_avg = mean(sensor_data(sensor_mask_index).Iz, ndims(sensor_data(sensor_mask_index).Iz));
        end
    end

end
    
% remove the non staggered particle velocity variables if not required
if ~flags.record_u_non_staggered
    switch kgrid.dim
        case 1
            sensor_data = rmfield(sensor_data, {'ux_non_staggered'});
        case 2
            sensor_data = rmfield(sensor_data, {'ux_non_staggered', 'uy_non_staggered'});
        case 3
            sensor_data = rmfield(sensor_data, {'ux_non_staggered', 'uy_non_staggered', 'uz_non_staggered'});
    end
end

% remove the time varying intensity if not required
if ~flags.record_I
    switch kgrid.dim
        case 1
            sensor_data = rmfield(sensor_data, {'Ix'});
        case 2
            sensor_data = rmfield(sensor_data, {'Ix', 'Iy'});
        case 3
            sensor_data = rmfield(sensor_data, {'Ix', 'Iy', 'Iz'});
    end
end

% remove the time varying pressure if not required
if ~flags.record_p
    sensor_data = rmfield(sensor_data, {'p'});
end

% remove unused variables
clear output_size ux_sgt uy_sgt uz_sgt sensor_mask_index;

%     % plot an example of the shifted velocity
%     figure;
%     plot((1:length(sensor_data.ux))*kgrid.dt, sensor_data.ux(1, :), 'k-s'); 
%     hold on;
%     plot((1:length(sensor_data.ux))*kgrid.dt + kgrid.dt/2, ux_sgt(1, :), 'r-s');
%     legend('Original', 'Shifted');