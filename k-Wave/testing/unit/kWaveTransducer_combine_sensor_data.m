function test_pass = kWaveTransducer_combine_sensor_data(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the output of the C++ and MATLAB codes when
%     using a k-Wave transducer as a sensor.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 13th January 2019
%     last update - 13th January 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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

% set comparison threshold
comparison_threshold = 2e-5;

% set pass variable
test_pass = true;

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 20;            % [grid points]
pml_y_size = 10;            % [grid points]
pml_z_size = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*pml_x_size;    % [grid points]
Ny = 128 - 2*pml_y_size;    % [grid points]
Nz = 64  - 2*pml_z_size;    % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]
medium.density = 1000;          % [kg/m^3]

% create the time array
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE SOURCE
% =========================================================================

% create source mask
source.p_mask = makeBall(Nx, Ny, Nz, round(25e-3/dx), Ny/2, Nz/2, 3)...
    + makeBall(Nx, Ny, Nz, round(8e-3/dx), Ny/4, Nz/2, 3);

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
source.p = source_strength .* toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 72;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points/voxels]
transducer.element_length = 12;     % length of each element [grid points/voxels]
transducer.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 25e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(21:52) = 1;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {...
    'DisplayMask', transducer.active_elements_mask, ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', 'single', ...
    'PlotScale', [-1/4, 1/4] * source_strength, ...
    'PlotSim', plot_simulations, ...
    };

% run the simulation using the transducer as sensor
sensor_data_trans = kspaceFirstOrder3D(kgrid, medium, source, transducer, input_args{:});    

% create a normal sensor mask
sensor.mask = double(transducer.active_elements_mask);

% run the simulation using a regular sensor
sensor_data_norm = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}); 

% combine the sensor data
sensor_data_norm = transducer.combine_sensor_data(sensor_data_norm);

% extract a single scan line from the sensor data using the current
% beamforming settings
scan_line_trans = transducer.scan_line(sensor_data_trans);
scan_line_norm  = transducer.scan_line(sensor_data_norm);

% trim the last few samples
trim_num = 20;
scan_line_trans   = scan_line_trans  (1:end - trim_num);
scan_line_norm    = scan_line_norm   (1:end - trim_num);
sensor_data_trans = sensor_data_trans(:, 1:end - trim_num);
sensor_data_norm  = sensor_data_norm (:, 1:end - trim_num);

% =========================================================================
% ERROR
% =========================================================================

% compute error
err = max(abs(sensor_data_trans(:) - sensor_data_norm(:))) / max(abs(sensor_data_trans(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons

    % plot the sensor data
    figure;
    subplot(3, 1, 1);
    imagesc(sensor_data_trans);
    colorbar;
    title('Transducer Sensor');
    
    subplot(3, 1, 2);
    imagesc(sensor_data_norm);
    colorbar;
    title('Regular Sensor');
    
    subplot(3, 1, 3)
    imagesc(100 * abs(sensor_data_trans - sensor_data_norm) ./ max(abs(sensor_data_trans(:))));
    colorbar;
    title('Error (%)');
    
    % plot the scan line
    figure;
    subplot(2, 1, 1);
    plot(kgrid.t_array(1:end - trim_num) * 1e6, scan_line_trans * 1e-6);
    hold on;
    plot(kgrid.t_array(1:end - trim_num) * 1e6, scan_line_norm * 1e-6, '--');
    xlabel('Time [\mus]');
    ylabel('Pressure [MPa]');
    title('Scan Line After Beamforming');

    subplot(2, 1, 2);
    plot(kgrid.t_array(1:end - trim_num) * 1e6, 100 * abs(scan_line_trans - scan_line_norm) ./ max(scan_line_trans(:)));
    xlabel('Time [\mus]');
    ylabel('Error (%)');
    
end