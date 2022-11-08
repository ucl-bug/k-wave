% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to scale
%     source terms to the correct units.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 15th February 2012
%     last update - 13th December 2018
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

% get the dimension size
N = kgrid.dim;

% check for non-uniform grid and give error for source terms that haven't
% yet been implemented
if (flags.nonuniform_grid && ( flags.source_uy || flags.source_uz || flags.transducer_source))
    disp('WARNING: source scaling not implemented for non-uniform grids with given source condition');
    return
end

% =========================================================================
% PRESSURE SOURCES
% =========================================================================

% apply k-space source correction expressed as a function of w
if flags.source_p && flags.use_w_source_correction_p
    source.p = source.p .* cos(2 * pi * source.p_frequency_ref * kgrid.dt/2);
end

% scale the input pressure by 1/c0^2 (to convert to units of density), then
% by 1/N (to split the input across the split density field). If the
% pressure is injected as a mass source, also scale the pressure by 
% 2*dt*c0/dx to account for the time step and convert to units of 
% [kg/(m^3 s)] 
if flags.source_p 
    if strcmp(source.p_mode, 'dirichlet')
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous
            % sound speed 
            source.p = source.p ./ (N .* c0.^2);
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for p_index = 1:length(source.p(:, 1))        
                source.p(p_index, :) = source.p(p_index, :) ./ (N .* c0(p_source_pos_index(p_index)).^2);
            end
            
        end
    else
        if flags.nonuniform_grid
            
            % create empty matrix
            grid_point_sep = zeros(size(kgrid.x));
            
            % compute averaged grid point seperation map, the interior
            % points are calculated using the average distance to all
            % connected grid points (the edge values are not calculated
            % assuming there are no source points in the PML)
            switch kgrid.dim
                case 1
                    grid_point_sep(2:end - 1) = ...
                        kgrid.x_size .* (kgrid.xn(3:end, 1) - kgrid.xn(1:end - 2, 1)) / 2;
                case 2
                    grid_point_sep(2:end - 1, 2:end - 1) = ...
                        kgrid.x_size .* (kgrid.xn(3:end, 2:end - 1) - kgrid.xn(1:end - 2, 2:end - 1)) / 4 + ...
                        kgrid.y_size .* (kgrid.yn(2:end - 1, 3:end) - kgrid.yn(2:end - 1, 1:end - 2)) / 4;
                case 3
                    grid_point_sep(2:end - 1, 2:end - 1, 2:end - 1) = ...
                        kgrid.x_size .* (kgrid.xn(3:end, 2:end - 1, 2:end - 1) - kgrid.xn(1:end - 2, 2:end - 1, 2:end - 1)) / 6 + ...
                        kgrid.y_size .* (kgrid.yn(2:end - 1, 3:end, 2:end - 1) - kgrid.yn(2:end - 1, 1:end - 2, 2:end - 1)) / 6 + ...
                        kgrid.z_size .* (kgrid.zn(2:end - 1, 2:end - 1, 3:end) - kgrid.zn(2:end - 1, 2:end - 1, 1:end - 2)) / 6;
            end
            
            % compute and apply scale parameter
            for p_index = 1:length(source.p(:, 1))
                if numel(c0) == 1
                    
                    % compute the scale parameter based on the homogeneous sound speed
                    source.p(p_index, :) = source.p(p_index, :) .* (2 .* dt ./ (N .* c0 .* grid_point_sep(p_source_pos_index(p_index))));
                    
                else
                    
                    % compute the scale parameter based on the sound speed at that position
                    source.p(p_index, :) = source.p(p_index, :) .* (2 .* dt ./ (N .* c0(p_source_pos_index(p_index)) .* grid_point_sep(p_source_pos_index(p_index))));
                    
                end
            end
            
            % clear unused variables
            clear grid_point_sep;
            
        else
            if numel(c0) == 1
                
                % compute the scale parameter based on the homogeneous
                % sound speed 
                source.p = source.p .* (2 .* dt ./ (N .* c0 .* kgrid.dx));
                
            else
                
                % compute the scale parameter seperately for each source
                % position based on the sound speed at that position
                for p_index = 1:length(source.p(:, 1))        
                    source.p(p_index, :) = source.p(p_index, :) .* (2 .* dt ./ (N .* c0(p_source_pos_index(p_index)) .* kgrid.dx));
                end
                
            end
        end
    end
end

% =========================================================================
% STRESS SOURCES
% =========================================================================

% scale the stress source by 1/N to divide amoungst the split field
% components, and if source.s_mode is not set to 'dirichlet', also scale by 
% 2*dt*c0/dx to account for the time step and convert to units of 
% [kg/(m^3 s)] (note dx is used in all dimensions)
if flags.source_sxx
    if strcmp(source.s_mode, 'dirichlet') || flags.source_p0
        source.sxx = source.sxx ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxx = source.sxx .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxx(:, 1))        
                source.sxx(s_index, :) = source.sxx(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end
if flags.source_syy
    if strcmp(source.s_mode, 'dirichlet') || flags.source_p0
        source.syy = source.syy ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.syy = source.syy .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.syy(:, 1))        
                source.syy(s_index, :) = source.syy(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end
if flags.source_szz
    if strcmp(source.s_mode, 'dirichlet') || flags.source_p0
        source.szz = source.szz ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.szz = source.szz .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.szz(:, 1))        
                source.szz(s_index, :) = source.szz(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end
if flags.source_sxy
    if strcmp(source.s_mode, 'dirichlet')
        source.sxy = source.sxy ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxy = source.sxy .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxy(:, 1))        
                source.sxy(s_index, :) = source.sxy(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end
if flags.source_sxz
    if strcmp(source.s_mode, 'dirichlet')
        source.sxz = source.sxz ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxz = source.sxz .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxz(:, 1))        
                source.sxz(s_index, :) = source.sxz(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end
if flags.source_syz
    if strcmp(source.s_mode, 'dirichlet')
        source.syz = source.syz ./ N;
    else
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.syz = source.syz .* (2 .* dt .* c0 ./ (N .* kgrid.dx));
            
        else
            
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.syz(:, 1))        
                source.syz(s_index, :) = source.syz(s_index, :) .* (2 .* dt .* c0(s_source_pos_index(s_index)) ./ (N .* kgrid.dx));
            end
            
        end         
    end
end

% =========================================================================
% VELOCITY SOURCES
% =========================================================================

% apply k-space source correction expressed as a function of w
if flags.use_w_source_correction_u
    if flags.source_ux
        source.ux = source.ux .* cos(2 * pi * source.u_frequency_ref * kgrid.dt/2);
    end
    if flags.source_uy
        source.uy = source.uy .* cos(2 * pi * source.u_frequency_ref * kgrid.dt/2);
    end
    if flags.source_uz
        source.uz = source.uz .* cos(2 * pi * source.u_frequency_ref * kgrid.dt/2);
    end
end

% if source.u_mode is not set to 'dirichlet', scale the x-direction
% velocity source terms by 2*dt*c0/dx to account for the time step and
% convert to units of [m/s^2] 
if flags.source_ux && ~strcmp(source.u_mode, 'dirichlet')
    
    if flags.nonuniform_grid
    
        % create empty matrix
        grid_point_sep = zeros(size(kgrid.x));
        
        % compute averaged grid point seperation map, the interior
        % points are calculated using the average distance to all
        % connected grid points (the edge values are not calculated
        % assuming there are no source points in the PML)
        grid_point_sep(2:end - 1, :, :) = kgrid.x_size .* (kgrid.xn(3:end, :, :) - kgrid.xn(1:end - 2, :, :)) / 2;
        
        % compute and apply scale parameter
        for u_index = 1:length(source.ux(:, 1))
            if numel(c0) == 1
                % compute the scale parameter based on the homogeneous sound speed
                source.ux(u_index, :) = source.ux(u_index, :) .* (2 .* c0 .* dt ./ (grid_point_sep(u_source_pos_index(u_index))));
            else
                % compute the scale parameter based on the sound speed at that position
                source.ux(u_index, :) = source.ux(u_index, :) .* (2 .* c0(u_source_pos_index(u_index)) .* dt ./ (grid_point_sep(u_source_pos_index(u_index))));
            end
        end
        
        % clear unused variables
        clear grid_point_sep;
    
    else
        
        if numel(c0) == 1
            
            % compute the scale parameter based on the homogeneous sound speed
            source.ux = source.ux .* (2 .* c0 .* dt ./ kgrid.dx);
            
        else
            
            % compute the scale parameter seperately for each source position
            % based on the sound speed at that position
            for u_index = 1:length(source.ux(:, 1))
                source.ux(u_index, :) = source.ux(u_index, :) .* (2 .* c0(u_source_pos_index(u_index)) .* dt ./ kgrid.dx);
            end
            
        end
    end
end

% if source.u_mode is not set to 'dirichlet', scale the y-direction
% velocity source terms by 2*dt*c0/dy to account for the time step and
% convert to units of [m/s^2] 
if flags.source_uy && ~strcmp(source.u_mode, 'dirichlet')
    if numel(c0) == 1
        
        % compute the scale parameter based on the homogeneous sound speed
        source.uy = source.uy .* (2 .* c0 .* dt ./ kgrid.dy);
        
    else
        
        % compute the scale parameter seperately for each source position
        % based on the sound speed at that position
        for u_index = 1:length(source.uy(:, 1))
            source.uy(u_index, :) = source.uy(u_index, :) .* (2 .* c0(u_source_pos_index(u_index)) .* dt ./ kgrid.dy);
        end
        
    end 
end 

% if source.u_mode is not set to 'dirichlet', scale the z-direction
% velocity source terms by 2*dt*c0/dz to account for the time step and
% convert to units of [m/s^2]  
if flags.source_uz && ~strcmp(source.u_mode, 'dirichlet') 
    if numel(c0) == 1
        
        % compute the scale parameter based on the homogeneous sound speed
        source.uz = source.uz .* (2 .* c0 .* dt ./ kgrid.dz);
        
    else
        
        % compute the scale parameter seperately for each source position
        % based on the sound speed at that position
        for u_index = 1:length(source.uz(:, 1))        
            source.uz(u_index, :) = source.uz(u_index, :) .* (2 .* c0(u_source_pos_index(u_index)) .* dt ./ kgrid.dz);
        end
        
    end
end

% =========================================================================
% TRANSDUCER SOURCE
% =========================================================================

% scale the transducer source term by 2*dt*c0/dx to account for the time
% step and convert to units of [m/s^2] 
if flags.transducer_source   
    if numel(c0) == 1
        transducer_input_signal = transducer_input_signal .* (2 .* c0 .* dt ./ kgrid.dx);
    else
        % compute the scale parameter based on the average sound speed at the
        % transducer positions (only one input signal is used to drive the
        % transducer)
        transducer_input_signal = transducer_input_signal .* (2 .* (mean(c0(u_source_pos_index))) .* dt ./ kgrid.dx);
    end
end

% clear subscript variables
clear N;