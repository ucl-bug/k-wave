% Equivalent Source Holography Example
%
% This examples demonstrates how to project a plane of measured data using
% kspaceFirstOrder3D by calculating an equivalent mass source (i.e.,
% source.p with source.p_mode = ’additive’) positioned in a parallel plane
% that recreates the measured data. It builds on the Simulations In Three
% Dimensions and Holographic Projections Using The Angular Spectrum Method
% examples.
%
% For a more detailed discussion of this example and the underlying
% techniques, see B. Treeby, F. Lucka, E. Martin, E. and B. Cox,
% "Equivalent-Source Acoustic Holography for Projecting Measured Ultrasound
% Fields Through Complex Media," IEEE Transactions on Ultrasonics,
% Ferroelectrics, and Frequency Control, vol. 65, no, 10, pp.1857-1864,
% 2018.
%
% author: Bradley Treeby
% date: 20th February 2019
% last update: 21st January 2020
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019-2020 Bradley Treeby

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

% =========================================================================
% DEFINE SIMULATION SETTINGS
% =========================================================================

% define grid properties
Nx = 44;                    % number of grid points in the x direction
Ny = 44;                    % number of grid points in the y direction
Nz = 44;                    % number of grid points in the z direction
dx = 1e-4;                  % grid point spacing [m]
c0 = 1500;                  % sound speed [m/s]

% define source properties
source_freq = 2e6;          % [Hz]
rect_x_radius = 8;          % [grid points]
rect_y_radius = 12;         % [grid points]
input_plane_z_index = 8;  	% [grid points]

% define simulation properties
cfl = 0.3;
t_end = 10e-6;               % [s]

% settings for calculating the equivalent source
grid_expansion = 0;
number_optimisation_steps = 20;

% =========================================================================
% SIMULATE MEASURED DATA
% =========================================================================

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% assign medium properties
medium.sound_speed = c0;

% compute sampling rates, forcing points-per-period to be an integer
points_per_wavelength = c0 / (source_freq * dx);
points_per_period = round(points_per_wavelength / cfl);

% compute corresponding time spacing
dt = 1 / (points_per_period * source_freq);    

% create the time array
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% define rectangular source mask
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(Nx/2 - rect_x_radius + 1:Nx/2 + rect_x_radius, Ny/2 - rect_y_radius + 1:Ny/2 + rect_y_radius, 1) = 1;

% define source signal as a continuous wave sinusoid
source.p = createCWSignals(kgrid.t_array, source_freq, 1, 0);

% define two planar sensor masks using opposing corners of a cuboid
sensor.mask = [1, 1,    input_plane_z_index, Nx, Ny,   input_plane_z_index; ...
               1, Ny/2, input_plane_z_index, Nx, Ny/2, Nz].';

% set the start time to only record the last three periods
sensor.record_start_index = kgrid.Nt - 3 * points_per_period + 1;
           
% assign input arguments
input_args = {...
    'PMLInside', false', ...
    'PMLSize', 'auto', ..., Nx, Ny 
    'DataCast', 'single', ...
    'DataRecast', true, ...
    };

% run k-Wave simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    
% extract the amplitude and phase from the time series data
[input_plane_amp, input_plane_phase] = extractAmpPhase(squeeze(sensor_data(1).p), 1/kgrid.dt, source_freq, ...
    'Dim', 3, ...
    'FFTPadding', 1, ...
    'Window', 'Rectangular');

[output_plane_amp, output_plane_phase] = extractAmpPhase(squeeze(sensor_data(2).p), 1/kgrid.dt, source_freq, ...
    'Dim', 3, ...
    'FFTPadding', 1, ...
    'Window', 'Rectangular');

% form data into matrix of complex values for use with angularSpectrumCW
input_plane_complex = input_plane_amp .* exp(1i * input_plane_phase);

% =========================================================================
% PROJECT MEASURED DATA USING A DIRICHLET BOUNDARY CONDITION
% =========================================================================

% assign source mask where data was measured
clear source;
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(:, :, input_plane_z_index) = 1;

% assign the measured data as a dirichlet boundary condition
source.p = createCWSignals(kgrid.t_array, source_freq, input_plane_amp(:), input_plane_phase(:));
source.p_mode = 'dirichlet';

% run k-Wave simulation to project measured data
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% extract amplitude from the time series data
proj_dirch = extractAmpPhase(squeeze(sensor_data(2).p), 1/kgrid.dt, source_freq, ...
    'Dim', 3, ...
    'FFTPadding', 1, ...
    'Window', 'Rectangular');

% =========================================================================
% PROJECT MEASURED DATA USING THE ANGULAR SPECTRUM METHOD
% =========================================================================

% run angular spectrum to project measured data
proj_asm = angularSpectrumCW(input_plane_complex, dx, ...
    (0:(Nz - input_plane_z_index)) * dx, source_freq, c0, ...
    'FFTLength', 512);

% extract amplitude from required plane
proj_asm = squeeze(abs(proj_asm(:, Ny/2, :)));

% =========================================================================
% CALCULATE EQUIVALENT SOURCE
% =========================================================================

% offset between the input plane and source plane [grid points]
source_offset = input_plane_z_index - 1;

% calculate equivalent source
[source_estimate, optim_params] = calculateMassSourceCW(input_plane_complex, dx, ...
    source_freq, c0, source_offset, grid_expansion, ...
    'NumSteps', number_optimisation_steps);

% =========================================================================
% PROJECT MEASURED DATA USING EQUIVALENT SOURCE
% =========================================================================

% assign source mask at the beginning of the grid
clear source;
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(:, :, 1) = 1;

% assign the measured data as a dirichlet boundary condition
source.p = createCWSignals(kgrid.t_array, source_freq, abs(source_estimate(:)), angle(source_estimate(:)));

% run k-Wave simulation to project measured data
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% extract amplitude from the time series data
proj_eqs = extractAmpPhase(squeeze(sensor_data(2).p), 1/kgrid.dt, source_freq, ...
    'Dim', 3, ...
    'FFTPadding', 1, ...
    'Window', 'Rectangular');

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the input data
figure;
subplot(1, 2, 1);
imagesc(input_plane_amp);
axis image;
colorbar;
title('Input Plane - Amplitude');

subplot(1, 2, 2);
imagesc(input_plane_phase);
axis image;
colorbar;
title('Input Plane - Phase');

% plot the reconstructed source plane
figure;
subplot(1, 2, 1);
imagesc(abs(source_estimate));
axis image;
colorbar;
title('Equivalent Source - Amplitude');

subplot(1, 2, 2);
imagesc(angle(source_estimate));
axis image;
colorbar;
title('Equivalent Source - Phase');

% plot the output data
figure;
subplot(2, 2, 1);
imagesc(output_plane_amp.');
axis image;
colorbar;
title('Reference Field');

subplot(2, 2, 2);
imagesc(100 * abs(output_plane_amp.' - proj_dirch.') ./ max(output_plane_amp(:)));
axis image;
colorbar;
title('Error [%] - Dirichlet Source');

subplot(2, 2, 3);
imagesc(100 * abs(output_plane_amp.' - proj_eqs.') ./ max(output_plane_amp(:)));
axis image;
colorbar;
title('Error [%] - Equivalent Source');

subplot(2, 2, 4);
imagesc(100 * abs(output_plane_amp.' - proj_asm.') ./ max(output_plane_amp(:)));
axis image;
colorbar;
title('Error [%] - Angular Spectrum');