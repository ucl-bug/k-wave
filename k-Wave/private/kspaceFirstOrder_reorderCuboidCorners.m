% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to reassign
%     the sensor data belonging to each set of cuboid corners from the
%     indexed sensor mask data. 
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

%#ok<*SAGROW>

% update command line status
disp('  reordering cuboid corners data...');

% set cuboid index variable
cuboid_start_pos = 1;

% loop through cuboid corners and for each recorded variable, reshape to
% [X, Y, Z, Y] or [X, Y, Z] instead of [sensor_index, T] or [sensor_index]
for cuboid_index = 1:size(record.cuboid_corners_list, 2)
    
    % set number of time points
    if flags.stream_to_disk
        cuboid_num_time_points = num_stream_time_points;
    else
        cuboid_num_time_points = num_recorded_time_points;
    end
    
    % get size of cuboid
    switch kgrid.dim
        case 1
            cuboid_size_x  = [1 + record.cuboid_corners_list(2, cuboid_index) - record.cuboid_corners_list(1, cuboid_index), 1];
            cuboid_size_xt = [cuboid_size_x(1), cuboid_num_time_points];
        case 2
            cuboid_size_x = [1 + record.cuboid_corners_list(3, cuboid_index) - record.cuboid_corners_list(1, cuboid_index),...
                             1 + record.cuboid_corners_list(4, cuboid_index) - record.cuboid_corners_list(2, cuboid_index)];
            cuboid_size_xt = [cuboid_size_x, cuboid_num_time_points];
        case 3
            cuboid_size_x = [1 + record.cuboid_corners_list(4, cuboid_index) - record.cuboid_corners_list(1, cuboid_index),...
                             1 + record.cuboid_corners_list(5, cuboid_index) - record.cuboid_corners_list(2, cuboid_index),...
                             1 + record.cuboid_corners_list(6, cuboid_index) - record.cuboid_corners_list(3, cuboid_index)];
            cuboid_size_xt = [cuboid_size_x, cuboid_num_time_points];
    end

    % set index and size variables
    cuboid_num_points = prod(cuboid_size_x);    
        
    if flags.record_p
        sensor_data_temp(cuboid_index).p = reshape(sensor_data.p(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt); 
    end

    if flags.record_p_max
        sensor_data_temp(cuboid_index).p_max = reshape(sensor_data.p_max(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
    end

    if flags.record_p_min
        sensor_data_temp(cuboid_index).p_min = reshape(sensor_data.p_min(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
    end

    if flags.record_p_rms
        sensor_data_temp(cuboid_index).p_rms = reshape(sensor_data.p_rms(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
    end

    if flags.record_u

        % x-dimension
        sensor_data_temp(cuboid_index).ux = reshape(sensor_data.ux(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        
        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).uy = reshape(sensor_data.uy(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).uz = reshape(sensor_data.uz(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end
        
    end

    if flags.record_u_non_staggered

        % x-dimension
        sensor_data_temp(cuboid_index).ux_non_staggered = reshape(sensor_data.ux_non_staggered(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        
        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).uy_non_staggered = reshape(sensor_data.uy_non_staggered(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).uz_non_staggered = reshape(sensor_data.uz_non_staggered(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end
        
    end    
    
    if flags.record_u_max

        % x-dimension
        sensor_data_temp(cuboid_index).ux_max = reshape(sensor_data.ux_max(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);

        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).uy_max = reshape(sensor_data.uy_max(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end    

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).uz_max = reshape(sensor_data.uz_max(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end
        
    end

    if flags.record_u_min

        % x-dimension
        sensor_data_temp(cuboid_index).ux_min = reshape(sensor_data.ux_min(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);

        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).uy_min = reshape(sensor_data.uy_min(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end    

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).uz_min = reshape(sensor_data.uz_min(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end
        
    end

    if flags.record_u_rms

        % x-dimension
        sensor_data_temp(cuboid_index).ux_rms = reshape(sensor_data.ux_rms(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);

        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).uy_rms = reshape(sensor_data.uy_rms(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end    

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).uz_rms = reshape(sensor_data.uz_rms(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end
        
    end

    if flags.record_I

        % x-dimension
        sensor_data_temp(cuboid_index).Ix = reshape(sensor_data.Ix(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        
        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).Iy = reshape(sensor_data.Iy(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).Iz = reshape(sensor_data.Iz(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1, :), cuboid_size_xt);
        end
        
    end
    
    if flags.record_I_avg

        % x-dimension
        sensor_data_temp(cuboid_index).Ix_avg = reshape(sensor_data.Ix_avg(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);

        % y-dimension if 2D or 3D
        if kgrid.dim > 1
            sensor_data_temp(cuboid_index).Iy_avg = reshape(sensor_data.Iy_avg(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end    

        % z-dimension if 3D
        if kgrid.dim > 2
            sensor_data_temp(cuboid_index).Iz_avg = reshape(sensor_data.Iz_avg(cuboid_start_pos:cuboid_start_pos + cuboid_num_points - 1), cuboid_size_x);
        end
        
    end
    
    % update cuboid index variable
    cuboid_start_pos = cuboid_start_pos + cuboid_num_points;
    
end

% assign max and final variables
if flags.record_p_final
    sensor_data_temp(1).p_final = sensor_data.p_final;
end

if flags.record_p_max_all
    sensor_data_temp(1).p_max_all = sensor_data.p_max_all;
end

if flags.record_p_min_all
    sensor_data_temp(1).p_min_all = sensor_data.p_min_all;
end

if flags.record_u_final

    % x-dimension
    sensor_data_temp(1).ux_final = sensor_data.ux_final;

    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data_temp(1).uy_final = sensor_data.uy_final;
    end

    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data_temp(1).uz_final = sensor_data.uz_final;
    end
    
end

if flags.record_u_max_all

    % x-dimension
    sensor_data_temp(1).ux_max_all = sensor_data.ux_max_all;

    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data_temp(1).uy_max_all = sensor_data.uy_max_all;
    end

    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data_temp(1).uz_max_all = sensor_data.uz_max_all;
    end
    
end

if flags.record_u_min_all

    % x-dimension
    sensor_data_temp(1).ux_min_all = sensor_data.ux_min_all;

    % y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data_temp(1).uy_min_all = sensor_data.uy_min_all;
    end

    % z-dimension if 3D
    if kgrid.dim > 2
        sensor_data_temp(1).uz_min_all = sensor_data.uz_min_all;
    end
    
end

% assign new sensor data to old
sensor_data = sensor_data_temp;

% clear any unused variables
clear cuboid_index cuboid_size_x cuboid_size_xt cuboid_start_pos cuboid_num_points cuboid_num_time_points;