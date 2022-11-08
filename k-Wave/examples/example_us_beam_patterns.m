% Simulating Ultrasound Beam Patterns Example
%
% This example shows how the nonlinear beam pattern from an ultrasound
% transducer can be modelled. It builds on the Defining An Ultrasound
% Transducer and Simulating Transducer Field Patterns examples. 
%
% author: Bradley Treeby
% date: 27th July 2011
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% simulation settings
DATA_CAST = 'single';       % set to 'single' or 'gpuArray-single' to speed up computations
MASK_PLANE = 'xy';          % set to 'xy' or 'xz' to generate the beam pattern in different planes
USE_STATISTICS = true;      % set to true to compute the rms or peak beam patterns, set to false to compute the harmonic beam patterns

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*PML_X_SIZE;    % [grid points]
Ny = 64 - 2*PML_Y_SIZE;     % [grid points]
Nz = 64 - 2*PML_Z_SIZE;     % [grid points]

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
medium.sound_speed = 1540;      % [m/s]
medium.density = 1000;          % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = 45e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;    	% [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 32;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points]
transducer.element_length = 12;     % length of each element [grid points]
transducer.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
switch MASK_PLANE
    case 'xy'
        
        % define mask
        sensor.mask(:, :, Nz/2) = 1;
        
        % store y axis properties        
        Nj = Ny;
        j_vec = kgrid.y_vec;
        j_label = 'y';
        
    case 'xz'
        
        % define mask
        sensor.mask(:, Ny/2, :) = 1;
        
        % store z axis properties
        Nj = Nz;
        j_vec = kgrid.z_vec;
        j_label = 'z';
        
end 

% set the record mode such that only the rms and peak values are stored
if USE_STATISTICS
    sensor.record = {'p_rms', 'p_max'};
end

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', transducer.all_elements_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotScale', [-1/2, 1/2] * source_strength};

% stream the data to disk in blocks of 100 if storing the complete time
% history 
if ~USE_STATISTICS
    input_args = [input_args {'StreamToDisk', 100}];
end

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});

% =========================================================================
% COMPUTE THE BEAM PATTERN USING SIMULATION STATISTICS
% =========================================================================

if USE_STATISTICS
    
    % reshape the returned rms and max fields to their original position
    sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nj]);
    sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);
    
    % plot the beam pattern using the pressure maximum
    figure;
    imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_max * 1e-6);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Total Beam Pattern Using Maximum Of Recorded Pressure');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;
    
    % plot the beam pattern using the pressure rms
    figure;
    imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, sensor_data.p_rms * 1e-6);
    xlabel([j_label '-position [mm]']);
    ylabel('x-position [mm]');
    title('Total Beam Pattern Using RMS Of Recorded Pressure');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;

    % end the example
    return
    
end

% =========================================================================
% COMPUTE THE BEAM PATTERN FROM THE AMPLITUDE SPECTRUM
% =========================================================================

% reshape the sensor data to its original position so that it can be
% indexed as sensor_data(x, j, t)
sensor_data = reshape(sensor_data, [Nx, Nj, kgrid.Nt]);

% compute the amplitude spectrum
[freq, amp_spect] = spect(sensor_data, 1/kgrid.dt, 'Dim', 3);

% compute the index at which the source frequency and its harmonics occur
[f1_value, f1_index] = findClosest(freq, tone_burst_freq);
[f2_value, f2_index] = findClosest(freq, 2 * tone_burst_freq);

% extract the amplitude at the source frequency and store
beam_pattern_f1 = amp_spect(:, :, f1_index);

% extract the amplitude at the second harmonic and store
beam_pattern_f2 = amp_spect(:, :, f2_index);       

% extract the integral of the total amplitude spectrum
beam_pattern_total = sum(amp_spect, 3);

% plot the beam patterns
figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, beam_pattern_f1 * 1e-3);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Beam Pattern At Source Fundamental');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [kPa]');
axis image;

figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, beam_pattern_f2 * 1e-3);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Beam Pattern At Second Harmonic');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [kPa]');
axis image;

figure;
imagesc(j_vec * 1e3, (kgrid.x_vec - min(kgrid.x_vec(:))) * 1e3, beam_pattern_total * 1e-6);
xlabel([j_label '-position [mm]']);
ylabel('x-position [mm]');
title('Total Beam Pattern Using Integral Of Recorded Pressure');
colormap(jet(256));
c = colorbar;
ylabel(c, 'Pressure [MPa]'); 
axis image;

% =========================================================================
% PLOT DIRECTIVITY PATTERN AT FOCUS
% =========================================================================

% compute the directivity at each of the harmonics
directivity_f1 = squeeze(beam_pattern_f1(round(transducer.focus_distance/dx), :));
directivity_f2 = squeeze(beam_pattern_f2(round(transducer.focus_distance/dx), :));

% normalise
directivity_f1 = directivity_f1 ./ max(directivity_f1(:));
directivity_f2 = directivity_f2 ./ max(directivity_f2(:));

% compute relative angles from transducer
if strcmp(MASK_PLANE, 'xy')
    horz_axis = ((1:Ny) - Ny/2) * dy;
else
    horz_axis = ((1:Nz) - Nz/2) * dz;
end
angles = 180 * atan2(horz_axis, transducer.focus_distance) / pi;

% plot the directivity
figure;
plot(angles, directivity_f1, 'k-', angles, directivity_f2, 'k--');
axis tight;
xlabel('Angle [deg]');
ylabel('Normalised Amplitude');
legend('Fundamental', 'Second Harmonic', 'Location', 'NorthWest'); 