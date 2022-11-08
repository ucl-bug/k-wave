% Running C++ Simulations Example
%
% This example demonstrates how to use the C++ versions of
% kspaceFirstOrder3D. Before use, the appropriate C++ codes must be
% downloaded from http://www.k-wave.org/download.php and placed in the
% binaries folder of the toolbox.
%
% author: Bradley Treeby
% date: 2nd December 2013
% last update: 1st August 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2019 Bradley Treeby

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

% modify this parameter to run the different examples
example_number = 1;
% 1: Save the input data to disk
% 2: Reload the output data from disk
% 3: Run the C++ simulation from MATLAB
% 4: Run the C++ simulation on a CUDA-enabled GPU from MATLAB

% input and output filenames (these must have the .h5 extension)
input_filename  = 'example_input.h5';
output_filename = 'example_output.h5';

% pathname for the input and output files
pathname = tempdir;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 256;                   % number of grid points in the x direction
Ny = 128;                   % number of grid points in the y direction
Nz = 64;                    % number of grid points in the z direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
dz = 0.1e-3;                % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% set the size of the PML
pml_size = 10;              % [grid points]

% define a scattering ball
ball_radius = 20;           % [grid points]
ball_x      = Nx/2 + 40;    % [grid points]
ball_y      = Ny/2;         % [grid points]
ball_z      = Nz/2;         % [grid points]
ball        = makeBall(Nx, Ny, Nz, ball_x, ball_y, ball_z, ball_radius);

% define the properties of the propagation medium
medium.sound_speed            = 1500*ones(Nx, Ny, Nz);	 % [m/s]
medium.sound_speed(ball == 1) = 1800;                    % [m/s]
medium.density                = 1000*ones(Nx, Ny, Nz);   % [kg/m^3]
medium.density(ball == 1)     = 1200;                    % [kg/m^3]
medium.alpha_coeff            = 0.75;                    % [dB/(MHz^y cm)]
medium.alpha_power            = 1.5;

% create the time array
Nt = 1200;                  % number of time steps
dt = 15e-9;                 % time step [s]
kgrid.setTime(Nt, dt);

% define a square source element facing in the x-direction
source_y_size = 60;         % [grid points]
source_z_size = 30;         % [grid points]
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(1 + pml_size, Ny/2 - source_y_size/2:Ny/2 + source_y_size/2, Nz/2 - source_z_size/2:Nz/2 + source_z_size/2) = 1;

% define a time varying sinusoidal source
source_freq     = 2e6;      % [Hz]
source_strength = 0.5e6;    % [Pa]
source.p        = source_strength * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, :, Nz/2) = 1;

% set the input arguments
input_args = {'PMLSize', pml_size};

switch example_number
    case 1
        
        % save the input data to disk and then exit
        kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', [pathname input_filename]);
        
        % display the required syntax to run the C++ simulation
        disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
        disp('Then, use the syntax shown below to run the simulation:');
        if isunix
            disp(['./kspaceFirstOrder-OMP -i ' pathname input_filename ' -o ' pathname output_filename ' --p_final --p_max']);
        else
            disp(['kspaceFirstOrder-OMP.exe -i "' pathname input_filename '" -o "' pathname output_filename '" --p_final --p_max']);
        end
            
        % stop the script
        return
        
    case 2
        
        % load output data from the C++ simulation
        sensor_data.p_final = h5read([pathname output_filename], '/p_final');
        sensor_data.p_max   = h5read([pathname output_filename], '/p_max');
        
    case 3

        % define the field parameters to record
        sensor.record = {'p_final', 'p_max'};

        % run the C++ simulation using the wrapper function
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});

    case 4

        % define the field parameters to record
        sensor.record = {'p_final', 'p_max'};

        % run the C++ simulation using the wrapper function
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});        
        
end

% =========================================================================
% VISUALISATION
% =========================================================================

% take an x-y slice through the final pressure output (this is recorded
% over the entire grid)
sensor_data.p_final = squeeze(sensor_data.p_final(:, :, Nz/2));

% reshape the maximum pressure output (this is recorded at the grid points
% specified by the sensor mask)
sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Ny]);

% add a display mask
ball_outline = makeCircle(Nx, Ny, ball_x, ball_y, ball_radius);
sensor_data.p_max  (ball_outline == 1 | source.p_mask(:, :, Nz/2) == 1) = max(sensor_data.p_max(:));
sensor_data.p_final(ball_outline == 1 | source.p_mask(:, :, Nz/2) == 1) = max(sensor_data.p_final(:));

% remove the pml
sensor_data.p_max   = sensor_data.p_max  (1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);
sensor_data.p_final = sensor_data.p_final(1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);

% get a suitable plot scale
x_vec = kgrid.x_vec(1 + pml_size:end - pml_size);
y_vec = kgrid.y_vec(1 + pml_size:end - pml_size);
[x_sc, scale, prefix] = scaleSI(max([x_vec; y_vec]));

% plot the final pressure field in the x-y plane
figure;
subplot(1, 2, 1);
imagesc(y_vec * scale, x_vec * scale, sensor_data.p_final, [-1, 1]*source_strength);
colormap(getColorMap);
xlabel(['y [' prefix 'm]']);
ylabel(['x [' prefix 'm]']);
axis image;
title('Final Pressure Field');

% plot the maximum pressure field in the x-y plane
subplot(1, 2, 2);
imagesc(y_vec * scale, x_vec * scale, sensor_data.p_max, [-2, 2]*source_strength);
colormap(getColorMap);
xlabel(['y [' prefix 'm]']);
ylabel(['x [' prefix 'm]']);
axis image;
title('Maximum Pressure');