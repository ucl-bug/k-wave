% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to cast
%     output variables back to double precision.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 1st September 2012
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2018 Bradley Treeby

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
disp('  recasting variables to double...');

% set recast string to cast variables back to double (the parallel
% computing toolbox requires the gather command to be used)
if strncmp(data_cast, 'gpuArray', 8)
    recast_pre = 'double(gather(';
    recast_post = '));';
else
    recast_pre = 'double(';
    recast_post = ');';
end

% time history of the acoustic pressure
if flags.record_p
    eval(['sensor_data.p = ' recast_pre 'sensor_data.p' recast_post]);
end

% maximum pressure
if flags.record_p_max
    eval(['sensor_data.p_max = ' recast_pre 'sensor_data.p_max' recast_post]);
end 

% minimum pressure
if flags.record_p_min
    eval(['sensor_data.p_min = ' recast_pre 'sensor_data.p_min' recast_post]);
end 

% rms pressure
if flags.record_p_rms
    eval(['sensor_data.p_rms = ' recast_pre 'sensor_data.p_rms' recast_post]);
end

% final acoustic pressure over all grid points
if flags.record_p_final
    eval(['sensor_data.p_final = ' recast_pre 'sensor_data.p_final' recast_post]);
end

% maximum pressure over all grid points
if flags.record_p_max_all
    eval(['sensor_data.p_max_all = ' recast_pre 'sensor_data.p_max_all' recast_post]);
end 

% minimum pressure over all grid points
if flags.record_p_min_all
    eval(['sensor_data.p_min_all = ' recast_pre 'sensor_data.p_min_all' recast_post]);
end 

% time history of the acoustic particle velocity
if flags.record_u
    eval(['sensor_data.ux = ' recast_pre 'sensor_data.ux' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy = ' recast_pre 'sensor_data.uy' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz = ' recast_pre 'sensor_data.uz' recast_post]);
    end
end

% time history of the acoustic particle velocity
if flags.record_u
    eval(['sensor_data.ux = ' recast_pre 'sensor_data.ux' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy = ' recast_pre 'sensor_data.uy' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz = ' recast_pre 'sensor_data.uz' recast_post]);
    end
end

% time history of the acoustic particle velocity on non-staggered grid
if flags.record_u_non_staggered
    eval(['sensor_data.ux_non_staggered = ' recast_pre 'sensor_data.ux_non_staggered' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_non_staggered = ' recast_pre 'sensor_data.uy_non_staggered' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_non_staggered = ' recast_pre 'sensor_data.uz_non_staggered' recast_post]);
    end
end

% time history of the split-field acoustic particle velocity on non-staggered grid
if flags.record_u_split_field
    eval(['sensor_data.ux_split_p = ' recast_pre 'sensor_data.ux_split_p' recast_post]);
    eval(['sensor_data.ux_split_s = ' recast_pre 'sensor_data.ux_split_s' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_split_p = ' recast_pre 'sensor_data.uy_split_p' recast_post]);
        eval(['sensor_data.uy_split_s = ' recast_pre 'sensor_data.uy_split_s' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_split_p = ' recast_pre 'sensor_data.uz_split_p' recast_post]);
        eval(['sensor_data.uz_split_s = ' recast_pre 'sensor_data.uz_split_s' recast_post]);
    end
end

% maximum particle velocity
if flags.record_u_max
    eval(['sensor_data.ux_max = ' recast_pre 'sensor_data.ux_max' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_max = ' recast_pre 'sensor_data.uy_max' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_max = ' recast_pre 'sensor_data.uz_max' recast_post]);
    end
end

% minimum particle velocity
if flags.record_u_min
    eval(['sensor_data.ux_min = ' recast_pre 'sensor_data.ux_min' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_min = ' recast_pre 'sensor_data.uy_min' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_min = ' recast_pre 'sensor_data.uz_min' recast_post]);
    end
end

% rms particle velocity
if flags.record_u_rms
    eval(['sensor_data.ux_rms = ' recast_pre 'sensor_data.ux_rms' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_rms = ' recast_pre 'sensor_data.uy_rms' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_rms = ' recast_pre 'sensor_data.uz_rms' recast_post]);
    end
end

% final particle velocity everywhere within medium
if flags.record_u_final
    eval(['sensor_data.ux_final = ' recast_pre 'sensor_data.ux_final' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_final = ' recast_pre 'sensor_data.uy_final' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_final = ' recast_pre 'sensor_data.uz_final' recast_post]);
    end
end

% maximum particle velocity over all grid points
if flags.record_u_max_all
    eval(['sensor_data.ux_max_all = ' recast_pre 'sensor_data.ux_max_all' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_max_all = ' recast_pre 'sensor_data.uy_max_all' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_max_all = ' recast_pre 'sensor_data.uz_max_all' recast_post]);
    end
end

% minimum particle velocity over all grid points
if flags.record_u_min_all
    eval(['sensor_data.ux_min_all = ' recast_pre 'sensor_data.ux_min_all' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_min_all = ' recast_pre 'sensor_data.uy_min_all' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_min_all = ' recast_pre 'sensor_data.uz_min_all' recast_post]);
    end
end

% object of the kWaveTransducer class is being used as a sensor
if flags.transducer_sensor
    eval(['sensor_data.transducer = ' recast_pre 'sensor_data.transducer' recast_post]);
end