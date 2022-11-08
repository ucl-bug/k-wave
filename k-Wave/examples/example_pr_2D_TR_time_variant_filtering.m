% Attenuation Compensation Using Time Variant Filtering Example
%
% This example demonstrates how the acoustic attenuation present in the
% photoacoustic forward problem can be compensated for using time-variant
% filtering. It builds on the 2D Time Reversal Reconstruction For A
% Circular Sensor and Attenuation Compensation Using Time Reversal
% examples.
%
% For a more detailed discussion of this example and the underlying
% techniques, see B. E. Treeby "Acoustic attenuation compensation in
% photoacoustic tomography using time-variant filtering," J. Biomed. Opt.,
% vol. 18, no. 3, p.036008, 2013.
%
% author: Bradley Treeby
% date: 14th August 2014
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2019 Bradley Treeby

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

clearvars;

% set to true to enable use of the parallel computing toolbox
USE_PARALLEL_COMPUTING_TOOLBOX = false;

% =========================================================================
% FORWARD SIMULATION
% =========================================================================

% define the size of the simulation grid and the PML
simulation_size = 512;                      % [grid points]
PML_size_forw = 20;                         % [grid points]
x = 52e-3;                                  % [m]
y = x;                                      % [m]

% reduce the number of grid points in Nx and Ny by the size of the PML so
% that using 'PMLInside' set to false will still give the correct
% simulation size 
Nx = simulation_size - 2 * PML_size_forw;   % [grid points]    
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium_non_absorbing.sound_speed = 1510;	% [m/s]

% create a duplicate of the propagation medium structure and append the
% absorption properties
medium = medium_non_absorbing;
medium.alpha_power = 1.5;      
medium.alpha_coeff = 3;                     % [dB/(MHz^y cm)]

% store maximum supported frequency
f_max = kgrid.k_max * medium.sound_speed / (2 * pi);

% load the shepp logan phantom (note if the simulation or PML sizes are
% changed, the loaded data will need to be resized)
load EXAMPLE_shepp_logan

% smooth the phantom and assign it to the initial pressure
shepp_logan = smooth(shepp_logan, true);
source.p0 = shepp_logan;

% define a circular Cartesian sensor mask
sensor_radius = 25e-3;                      % [m]
sensor_points = 200;
cart_sensor_mask = makeCartCircle(sensor_radius, sensor_points);
sensor.mask = cart_sensor_mask;

% create the time array used for the simulation, with t_max defined using
% Huygens' principle to avoid artifact trapping in the reconstruction
t_max = 2 * sensor_radius / medium.sound_speed;
kgrid.makeTime(medium.sound_speed, [], t_max);

% set the input options, switching off the smoothing (the input has already
% been smoothed), setting the PML to be outside the defined grid, casting
% to 'single' to speed up the example, and switching off visualisation
input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_forw, 'DataCast', 'single', 'PlotSim', true};

% run the forward simulation
sensor_data_lossy = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise to the recorded sensor data
signal_to_noise_ratio = 40;                 % [dB]
sensor_data_lossy = addNoise(sensor_data_lossy, signal_to_noise_ratio, 'peak');

% =========================================================================
% IMAGE RECONSTRUCTION WITHOUT ATTENUATION COMPENSATION
% =========================================================================

% create a second computation grid for the reconstruction to avoid the
% inverse crime
PML_size_recon = 25;                        % [grid points]
Nx = simulation_size - 2 * PML_size_recon;  % [grid points]
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]
kgrid_recon = kWaveGrid(Nx, dx, Ny, dy);

% use the same time array for the reconstruction
kgrid_recon.setTime(kgrid.Nt, kgrid.dt);

% remove the initial pressure field from the source structure
source = rmfield(source, 'p0');

% create a continuous binary sensor mask with the same radius as the
% Cartesian sensor mask used in the forward simulation
radius_gp = round(sensor_radius/kgrid_recon.dx);
binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, floor(kgrid_recon.Nx/2) + 1, floor(kgrid_recon.Ny/2) + 1, radius_gp);

% assign the sensor mask to the sensor structure
sensor.mask = binary_sensor_mask;

% interpolate the simulated sensor data onto a continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data_lossy, cart_sensor_mask, binary_sensor_mask, 'linear');

% re-assign the input options (the PML_size has changed)
input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_recon, 'DataCast', 'single', 'PlotSim', true};

% run the time-reversal reconstruction using the non absorbing medium
p0_recon = kspaceFirstOrder2D(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});

% =========================================================================
% IMAGE RECONSTRUCTION WITH ATTENUATION COMPENSATION - FIXED
% =========================================================================

% correct for acoustic attenuation using time variant filtering regularised
% by a Tukey window with a fixed cutoff frequency
sensor_data_comp = attenComp(sensor_data_lossy, kgrid.dt, medium.sound_speed, ...
    medium.alpha_coeff, medium.alpha_power, 'FilterCutoff', [3e6, 3e6]);

% interpolate the simulated sensor data onto a continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data_comp, cart_sensor_mask, binary_sensor_mask, 'linear');
            
% run the time-reversal reconstruction
p0_recon_comp_fixed = kspaceFirstOrder2D(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});

% =========================================================================
% IMAGE RECONSTRUCTION WITH ATTENUATION COMPENSATION - AVERAGE
% =========================================================================

% correct for acoustic attenuation using time variant filtering regularised
% by a Tukey window with a time-variant cutoff frequency based on the
% average time-frequency distribution of the signals
sensor_data_comp = attenComp(sensor_data_lossy, kgrid.dt, medium.sound_speed, ...
    medium.alpha_coeff, medium.alpha_power, 'FrequencyMultiplier', 3, 'Plot', true);

% interpolate the simulated sensor data onto a continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data_comp, cart_sensor_mask, binary_sensor_mask, 'linear');

% run the time-reversal reconstruction
p0_recon_comp_avg = kspaceFirstOrder2D(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});

% =========================================================================
% IMAGE RECONSTRUCTION WITH ATTENUATION COMPENSATION - INDIVIDUAL
% =========================================================================

% correct for acoustic attenuation using time variant filtering regularised
% by a Tukey window with a time-variant cutoff frequency based on the
% time-frequency distribution for each signal
sensor_data_comp = zeros(size(sensor_data_lossy));

if ~USE_PARALLEL_COMPUTING_TOOLBOX
    
    % update command line
    disp('Appying Time Variant Filter... ');
    
    % loop through signals
    for index = 1:size(sensor_data_lossy, 1)

        % correct signal
        sensor_data_comp(index, :) = attenComp(sensor_data_lossy(index, :), kgrid.dt,...
            medium.sound_speed, medium.alpha_coeff, medium.alpha_power, 'FrequencyMultiplier', 3);

        % display signals
        if index == 1
            figure;
        end
        plot(kgrid.t_array * 1e6, sensor_data_lossy(index, :), 'k-', ...
             kgrid.t_array * 1e6, sensor_data_comp(index, :), 'r-');
        title(['Signal ' num2str(index) ' of ' num2str(sensor_points)]);
        drawnow;

    end
    
    % update command line
    disp(['  completed in ' scaleTime(toc)]);

else

    % update command line
    disp('Appying Time Variant Filter... ');    
    
    % start matlab pool
    parpool;
    
    % loop through signals using parfor
    parfor index = 1:size(sensor_data_lossy, 1)
        
        % correct signal
        sensor_data_comp(index, :) = attenComp(sensor_data_lossy(index, :), kgrid.dt,...
            medium.sound_speed, medium.alpha_coeff, medium.alpha_power, ...
            'FrequencyMultiplier', 3, 'DisplayUpdates', false); %#ok<PFBNS>
        
    end
    
    % delete matlab pool
    delete(gcp('nocreate'));
    
    % update command line
    disp(['  completed in ' scaleTime(toc)]);
    
end

% interpolate the simulated sensor data onto a continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data_comp, cart_sensor_mask, binary_sensor_mask, 'linear');

% run the time-reversal reconstruction
p0_recon_comp_tvf = kspaceFirstOrder2D(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% set the pressure outside the sensor mask to be zero
binary_sensor_map = makeDisc(kgrid_recon.Nx, kgrid_recon.Ny, kgrid_recon.Nx/2, kgrid_recon.Ny/2, radius_gp - 50);
p0_recon(binary_sensor_map ~= 1) = 0;
p0_recon_comp_tvf(binary_sensor_map ~= 1) = 0;
p0_recon_comp_avg(binary_sensor_map ~= 1) = 0;
p0_recon_comp_fixed(binary_sensor_map ~= 1) = 0;

% set the plot limits
x_lim = [-20, 20];
y_lim = [-20, 20];

% ----

% plot the original shepp logan phantom
figure;
subplot(1, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, shepp_logan, [0, 1]);
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Initial Pressure Distribution');
scaleFig(1.2, 0.8);

% ----

% plot the reconstruction with no filter
figure;
subplot(1, 2, 1);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon, [0, 1]);
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Reconstruction - No Filter');

% plot a profile
subplot(1, 2, 2);
plot(kgrid.y_vec * 1e3, shepp_logan(end/2,:), 'k');
hold on;
plot(kgrid_recon.y_vec * 1e3, p0_recon(end/2,:), 'r-');
axis square;
set(gca, 'XLim', x_lim, 'YLim', [-0.1, 1.1]);
ylabel('Pressure Magnitude [au]');
xlabel('y-position [mm]');
title('Profile');
scaleFig(1.2, 0.8);
 
% ----

% plot the reconstruction with fixed tvf filter
figure;
subplot(1, 2, 1);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_comp_fixed, [0, 1]);
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Reconstruction - Fixed Time Variant Filter');

% plot a profile
subplot(1, 2, 2);
plot(kgrid.y_vec * 1e3, shepp_logan(end/2,:), 'k');
hold on;
plot(kgrid_recon.y_vec * 1e3, p0_recon_comp_fixed(end/2,:), 'r-');
axis square;
set(gca, 'XLim', x_lim, 'YLim', [-0.1, 1.1]);
ylabel('Pressure Magnitude [au]');
xlabel('y-position [mm]');
title('Profile');
scaleFig(1.2, 0.8);

% ----

% plot the reconstruction with average tvf filter
figure;
subplot(1, 2, 1);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_comp_avg, [0, 1]);
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Reconstruction - Average Time Variant Filter');

% plot a profile
subplot(1, 2, 2);
plot(kgrid.y_vec * 1e3, shepp_logan(end/2,:), 'k');
hold on;
plot(kgrid_recon.y_vec * 1e3, p0_recon_comp_avg(end/2,:), 'r-');
axis square;
set(gca, 'XLim', x_lim, 'YLim', [-0.1, 1.1]);
ylabel('Pressure Magnitude [au]');
xlabel('y-position [mm]');
title('Profile');
scaleFig(1.2, 0.8);

% ----

% plot the reconstruction with tvf filter
figure;
subplot(1, 2, 1);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_comp_tvf, [0, 1]);
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Reconstruction - Time Variant Filter');

% plot a profile
subplot(1, 2, 2);
plot(kgrid.y_vec * 1e3, shepp_logan(end/2,:), 'k');
hold on;
plot(kgrid_recon.y_vec * 1e3, p0_recon_comp_tvf(end/2,:), 'r-');
axis square;
set(gca, 'XLim', x_lim, 'YLim', [-0.1, 1.1]);
ylabel('Pressure Magnitude [au]');
xlabel('y-position [mm]');
title('Profile');
scaleFig(1.2, 0.8);