function test_pass = kspaceFirstOrder2D_compare_binary_and_cuboid_sensor_mask(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the simulation results using a binary and
%     cuboid sensor mask
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th July 2014
%     last update - 19th February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;  
end

% set pass variable
test_pass = true;

% set additional literals to give further permutations of the test
COMPARISON_THRESH = 1e-15;
PML_INSIDE        = true;    
DATA_CAST         = 'off';

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x direction
Ny = 128;           % number of grid points in the y direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny);	% [m/s]
medium.sound_speed(1:Nx/2, :) = 1800;    % [m/s]
medium.density = 1000*ones(Nx, Ny);     % [kg/m^3]
medium.density(:, Ny/4:end) = 1200;      % [kg/m^3]

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2;

% define list of cuboid corners using two intersecting cuboids
cuboid_corners = [40, 10, 90, 65;...
                  10, 60, 50, 70].';
sensor.mask = cuboid_corners;

% set the variables to record
sensor.record = {'p', 'p_max', 'p_min', 'p_rms', 'p_max_all', 'p_min_all', 'p_final', ...
                 'u', 'u_max', 'u_min', 'u_rms', 'u_max_all', 'u_min_all', 'u_final',...
                 'I', 'I_avg'};

% set input options
input_args = {'PMLInside', PML_INSIDE, 'DataCast', DATA_CAST, 'PlotSim', plot_simulations};
             
% run the simulation
sensor_data_cuboids = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% ------------------------

% create a binary mask for display from the list of corners
sensor.mask = false(size(kgrid.k));
cuboid_index = 1;
sensor.mask(cuboid_corners(1, cuboid_index):cuboid_corners(3, cuboid_index),...
            cuboid_corners(2, cuboid_index):cuboid_corners(4, cuboid_index)) = 1;

% run the simulation
sensor_data_comp1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
         
% compute the error from the first cuboid
L_inf_p      = max(abs(sensor_data_cuboids(cuboid_index).p(:)      - sensor_data_comp1.p(:)))     / max(abs(sensor_data_comp1.p(:)));
L_inf_p_max  = max(abs(sensor_data_cuboids(cuboid_index).p_max(:)  - sensor_data_comp1.p_max(:))) / max(abs(sensor_data_comp1.p_max(:)));
L_inf_p_min  = max(abs(sensor_data_cuboids(cuboid_index).p_min(:)  - sensor_data_comp1.p_min(:))) / max(abs(sensor_data_comp1.p_min(:)));
L_inf_p_rms  = max(abs(sensor_data_cuboids(cuboid_index).p_rms(:)  - sensor_data_comp1.p_rms(:))) / max(abs(sensor_data_comp1.p_rms(:)));

L_inf_ux      = max(abs(sensor_data_cuboids(cuboid_index).ux(:)      - sensor_data_comp1.ux(:)))     / max(abs(sensor_data_comp1.ux(:)));
L_inf_ux_max  = max(abs(sensor_data_cuboids(cuboid_index).ux_max(:)  - sensor_data_comp1.ux_max(:))) / max(abs(sensor_data_comp1.ux_max(:)));
L_inf_ux_min  = max(abs(sensor_data_cuboids(cuboid_index).ux_min(:)  - sensor_data_comp1.ux_min(:))) / max(abs(sensor_data_comp1.ux_min(:)));
L_inf_ux_rms  = max(abs(sensor_data_cuboids(cuboid_index).ux_rms(:)  - sensor_data_comp1.ux_rms(:))) / max(abs(sensor_data_comp1.ux_rms(:)));

L_inf_uy      = max(abs(sensor_data_cuboids(cuboid_index).uy(:)      - sensor_data_comp1.uy(:)))     / max(abs(sensor_data_comp1.uy(:)));
L_inf_uy_max  = max(abs(sensor_data_cuboids(cuboid_index).uy_max(:)  - sensor_data_comp1.uy_max(:))) / max(abs(sensor_data_comp1.uy_max(:)));
L_inf_uy_min  = max(abs(sensor_data_cuboids(cuboid_index).uy_min(:)  - sensor_data_comp1.uy_min(:))) / max(abs(sensor_data_comp1.uy_min(:)));
L_inf_uy_rms  = max(abs(sensor_data_cuboids(cuboid_index).uy_rms(:)  - sensor_data_comp1.uy_rms(:))) / max(abs(sensor_data_comp1.uy_rms(:)));

% compute the error from the total variables
L_inf_p_max_all  = max(abs(sensor_data_cuboids(cuboid_index).p_max_all(:)  - sensor_data_comp1.p_max_all(:)))  / max(abs(sensor_data_comp1.p_max_all(:)));
L_inf_ux_max_all = max(abs(sensor_data_cuboids(cuboid_index).ux_max_all(:) - sensor_data_comp1.ux_max_all(:))) / max(abs(sensor_data_comp1.ux_max_all(:)));
L_inf_uy_max_all = max(abs(sensor_data_cuboids(cuboid_index).uy_max_all(:) - sensor_data_comp1.uy_max_all(:))) / max(abs(sensor_data_comp1.uy_max_all(:)));

L_inf_p_min_all  = max(abs(sensor_data_cuboids(cuboid_index).p_min_all(:)  - sensor_data_comp1.p_min_all(:)))  / max(abs(sensor_data_comp1.p_min_all(:)));
L_inf_ux_min_all = max(abs(sensor_data_cuboids(cuboid_index).ux_min_all(:) - sensor_data_comp1.ux_min_all(:))) / max(abs(sensor_data_comp1.ux_min_all(:)));
L_inf_uy_min_all = max(abs(sensor_data_cuboids(cuboid_index).uy_min_all(:) - sensor_data_comp1.uy_min_all(:))) / max(abs(sensor_data_comp1.uy_min_all(:)));

L_inf_p_final  = max(abs(sensor_data_cuboids(cuboid_index).p_final(:)  - sensor_data_comp1.p_final(:)))  / max(abs(sensor_data_comp1.p_final(:)));
L_inf_ux_final = max(abs(sensor_data_cuboids(cuboid_index).ux_final(:) - sensor_data_comp1.ux_final(:))) / max(abs(sensor_data_comp1.ux_final(:)));
L_inf_uy_final = max(abs(sensor_data_cuboids(cuboid_index).uy_final(:) - sensor_data_comp1.uy_final(:))) / max(abs(sensor_data_comp1.uy_final(:)));

% get maximum error
L_inf_max = max([L_inf_p, L_inf_p_max, L_inf_p_min, L_inf_p_rms,... 
    L_inf_ux, L_inf_ux_max, L_inf_ux_min, L_inf_ux_rms, ...
    L_inf_uy, L_inf_uy_max, L_inf_uy_min, L_inf_uy_rms, ...
    L_inf_p_max_all, L_inf_ux_max_all, L_inf_uy_max_all, ...
    L_inf_p_min_all, L_inf_ux_min_all, L_inf_uy_min_all, ...
    L_inf_p_final, L_inf_ux_final, L_inf_uy_final]);

% compute pass
if (L_inf_max > COMPARISON_THRESH)
    test_pass = false;
end

% ------------------------

% create a binary mask for display from the list of corners
sensor.mask = false(size(kgrid.k));
cuboid_index = 2;
sensor.mask(cuboid_corners(1, cuboid_index):cuboid_corners(3, cuboid_index),...
            cuboid_corners(2, cuboid_index):cuboid_corners(4, cuboid_index)) = 1;

% run the simulation
sensor_data_comp2 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
         
% compute the error from the second cuboid
L_inf_p      = max(abs(sensor_data_cuboids(cuboid_index).p(:)      - sensor_data_comp2.p(:)))     / max(abs(sensor_data_comp2.p(:)));
L_inf_p_max  = max(abs(sensor_data_cuboids(cuboid_index).p_max(:)  - sensor_data_comp2.p_max(:))) / max(abs(sensor_data_comp2.p_max(:)));
L_inf_p_min  = max(abs(sensor_data_cuboids(cuboid_index).p_min(:)  - sensor_data_comp2.p_min(:))) / max(abs(sensor_data_comp2.p_min(:)));
L_inf_p_rms  = max(abs(sensor_data_cuboids(cuboid_index).p_rms(:)  - sensor_data_comp2.p_rms(:))) / max(abs(sensor_data_comp2.p_rms(:)));

L_inf_ux      = max(abs(sensor_data_cuboids(cuboid_index).ux(:)      - sensor_data_comp2.ux(:)))     / max(abs(sensor_data_comp2.ux(:)));
L_inf_ux_max  = max(abs(sensor_data_cuboids(cuboid_index).ux_max(:)  - sensor_data_comp2.ux_max(:))) / max(abs(sensor_data_comp2.ux_max(:)));
L_inf_ux_min  = max(abs(sensor_data_cuboids(cuboid_index).ux_min(:)  - sensor_data_comp2.ux_min(:))) / max(abs(sensor_data_comp2.ux_min(:)));
L_inf_ux_rms  = max(abs(sensor_data_cuboids(cuboid_index).ux_rms(:)  - sensor_data_comp2.ux_rms(:))) / max(abs(sensor_data_comp2.ux_rms(:)));

L_inf_uy      = max(abs(sensor_data_cuboids(cuboid_index).uy(:)      - sensor_data_comp2.uy(:)))     / max(abs(sensor_data_comp2.uy(:)));
L_inf_uy_max  = max(abs(sensor_data_cuboids(cuboid_index).uy_max(:)  - sensor_data_comp2.uy_max(:))) / max(abs(sensor_data_comp2.uy_max(:)));
L_inf_uy_min  = max(abs(sensor_data_cuboids(cuboid_index).uy_min(:)  - sensor_data_comp2.uy_min(:))) / max(abs(sensor_data_comp2.uy_min(:)));
L_inf_uy_rms  = max(abs(sensor_data_cuboids(cuboid_index).uy_rms(:)  - sensor_data_comp2.uy_rms(:))) / max(abs(sensor_data_comp2.uy_rms(:)));

% get maximum error
L_inf_max = max([L_inf_p, L_inf_p_max, L_inf_p_min, L_inf_p_rms,... 
    L_inf_ux, L_inf_ux_max, L_inf_ux_min, L_inf_ux_rms, ...
    L_inf_uy, L_inf_uy_max, L_inf_uy_min, L_inf_uy_rms]);

% compute pass
if (L_inf_max > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% PLOT COMPARISONS
% =========================================================================

if plot_comparisons

    % plot the simulated sensor data
    figure;
    subplot(3, 2, 1);
    imagesc(reshape(sensor_data_cuboids(1).p, [], size(sensor_data_comp1.p, 2)), [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 1');
    colorbar;
    
    subplot(3, 2, 2);
    imagesc(reshape(sensor_data_cuboids(2).p, [], size(sensor_data_comp1.p, 2)), [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 2');
    colorbar;
    
    subplot(3, 2, 3);
    imagesc(sensor_data_comp1.p, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 1 - Comparison');
    colorbar;
    
    subplot(3, 2, 4);
    imagesc(sensor_data_comp2.p, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 2 - Comparison');
    colorbar;    
    
    subplot(3, 2, 5);
    imagesc(reshape(sensor_data_cuboids(1).p, [], size(sensor_data_comp1.p, 2)) - sensor_data_comp1.p, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 1 - Difference');
    colorbar;
    
    subplot(3, 2, 6);
    imagesc(reshape(sensor_data_cuboids(2).p, [], size(sensor_data_comp1.p, 2)) - sensor_data_comp2.p, [-1, 1]);
    colormap(getColorMap);
    ylabel('Sensor Position');
    xlabel('Time Step');
    title('Cuboid 2 - Difference');
    colorbar;      
    
end