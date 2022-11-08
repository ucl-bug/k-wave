function test_pass = kspaceFirstOrder2D_check_nonstaggered_particle_velocity(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to check that the output of u_non_staggered with u using
%       kspaceFirstOrder2D
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 31st July 2013
%       last update - 21st August 2014
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

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

% set pass variable
test_pass = true;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.density = 1000;      % [kg/m^3]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5; 

% create time array
CFL = 0.1;
t_end = 6e-6;
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [au]
disc_x_pos = Nx/2;  % [grid points]
disc_y_pos = Ny/2;  % [grid points]
disc_radius = 5;    % [grid points]
source.p0 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define four sensor points centered about source.p0
sensor_radius = 40; % [grid points]
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2 + sensor_radius + 1, Ny/2) = 1;

% set the acoustic variables that are recorded
sensor.record = {'p', 'u', 'u_non_staggered'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% find the maximum in each trace
[max1, ind1] = max(sensor_data.ux_non_staggered(1, :));
[max2, ind2] = max(sensor_data.ux(1, :));
[max3, ind3] = max(sensor_data.ux_non_staggered(2, :));
[max4, ind4] = max(sensor_data.ux(2, :));

% make sure the maximum are increasing in time and decreasing in magnitude
if (max1 < max2) || (max2 < max3) || (max3 < max4) ...
        || (ind1 > ind2) || (ind2 > ind3) || (ind3 > ind4)
    test_pass = false;
end

% run the simulation again using a Cartesian sensor mask
sensor.mask = grid2cart(kgrid, sensor.mask);
sensor_data_cart = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', plot_simulations);

% find the maximum in each trace
[max1, ind1] = max(sensor_data.ux_non_staggered(1, :));
[max2, ind2] = max(sensor_data.ux(1, :));
[max3, ind3] = max(sensor_data.ux_non_staggered(2, :));
[max4, ind4] = max(sensor_data.ux(2, :));

% make sure the maximum are increasing in time and decreasing in magnitude
if (max1 < max2) || (max2 < max3) || (max3 < max4) ...
        || (ind1 > ind2) || (ind2 > ind3) || (ind3 > ind4)
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons

    % plot the simulated sensor data
    [t, t_sc, t_prefix] = scaleSI(kgrid.t_array(end));
    mx = 5e-7;

    % plot binary sensor data
    figure;
    plot(t_sc*kgrid.t_array, sensor_data.ux_non_staggered(1, :), 'g-');
    hold on;
    plot(t_sc*kgrid.t_array, sensor_data.ux(1, :), 'r-');
    plot(t_sc*kgrid.t_array, sensor_data.ux_non_staggered(2, :), 'k-');    
    plot(t_sc*kgrid.t_array, sensor_data.ux(2, :), 'b-');
    set(gca, 'YLim', [-mx, mx], 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('ux'); 
    legend('1', '2', '3', '4');
    title('Binary sensor mask')

    % plot cartesian sensor data
    figure;
    plot(t_sc*kgrid.t_array, sensor_data_cart.ux_non_staggered(1, :), 'g-');
    hold on;    
    plot(t_sc*kgrid.t_array, sensor_data_cart.ux(1, :), 'r-');
    plot(t_sc*kgrid.t_array, sensor_data_cart.ux_non_staggered(2, :), 'k-');    
    plot(t_sc*kgrid.t_array, sensor_data_cart.ux(2, :), 'b-');
    set(gca, 'YLim', [-mx, mx], 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('ux'); 
    legend('1', '2', '3', '4');
    title('Cartesian sensor mask');

    % plot difference for good measure
    figure;
    plot(t_sc*kgrid.t_array, abs(sensor_data_cart.ux_non_staggered(1, :) - sensor_data.ux_non_staggered(1, :)), 'g-');
    hold on;    
    plot(t_sc*kgrid.t_array, abs(sensor_data_cart.ux(1, :) - sensor_data.ux(1, :)), 'r-');
    plot(t_sc*kgrid.t_array, abs(sensor_data_cart.ux_non_staggered(2, :) - sensor_data.ux_non_staggered(2, :)), 'k-');    
    plot(t_sc*kgrid.t_array, abs(sensor_data_cart.ux(2, :) - sensor_data.ux(2, :)), 'b-');
    set(gca, 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('ux'); 
    legend('1', '2', '3', '4');
    title('Difference between binary and cartesian sensor mask');
    
end
