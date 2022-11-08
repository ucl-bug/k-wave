% Saving Input Files In Parts Example
%
% This example demonstrates how to save the HDF5 input files required by
% the C++ code in parts. It builds on the Running C++ Simulations Example. 
%
% author: Bradley Treeby
% date: 3rd December 2013
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
% 1: Save the input data to disk in parts
% 2: Reload the output data from disk

% input and output filenames (these must have the .h5 extension)
input_filename  = 'example_input.h5';
output_filename = 'example_output.h5';

% pathname for the input and output files
pathname = tempdir;

% remove input file if it already exists
if example_number == 1 && exist([pathname input_filename], 'file')
    delete([pathname input_filename]);
end

% load HDF5 constants
run([getkWavePath 'private/getH5Literals']);

% =========================================================================
% SIMULATION SETTINGS
% =========================================================================

% set the properties of the computational grid
Nx = 256;                   % number of grid points in the x direction
Ny = 128;                   % number of grid points in the y direction
Nz = 64;                    % number of grid points in the z direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
dz = 0.1e-3;                % grid point spacing in the z direction [m]
Nt = 1200;                  % number of time steps
dt = 15e-9;                 % time step [s]

% set the properties of the perfectly matched layer 
pml_x_size  = 10;           % [grid points]
pml_y_size  = 10;           % [grid points]
pml_z_size  = 10;           % [grid points]
pml_x_alpha = 2;            % [Nepers/grid point]
pml_y_alpha = 2;            % [Nepers/grid point]
pml_z_alpha = 2;            % [Nepers/grid point]

% define a scattering ball
ball_radius = 20;           % [grid points]
ball_x      = Nx/2 + 40;    % [grid points]
ball_y      = Ny/2;         % [grid points]
ball_z      = Nz/2;         % [grid points]

% define the properties of the medium
c0_background   = 1500;     % [kg/m^3]
c0_ball         = 1800;     % [kg/m^3]
rho0_background = 1000;     % [kg/m^3]
rho0_ball       = 1200;     % [kg/m^3]
alpha_coeff     = 0.75;     % [dB/(MHz^y cm)]
alpha_power     = 1.5;

% define a the properties of a single square source element facing in the
% x-direction 
source_y_size   = 60;       % [grid points]
source_z_size   = 30;       % [grid points]
source_freq     = 2e6;      % [Hz]
source_strength = 0.5e6;    % [Pa]

% =========================================================================
% WRITE THE INPUT FILE
% =========================================================================

if example_number == 1

    % ---------------------------------------------------------------------
    % WRITE THE MEDIUM PARAMETERS
    % ---------------------------------------------------------------------
    
    % update command line status
    tic; fprintf('Writing medium parameters... ');

    % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

    % create the scattering ball and density matrix
    ball       = makeBall(Nx, Ny, Nz, ball_x, ball_y, ball_z, ball_radius, [], true);
    rho0       = rho0_background * ones(Nx, Ny, Nz, 'single');   
    rho0(ball) = rho0_ball;

    % make sure the input is in the correct data format
    eval(['rho0 = ' MATRIX_DATA_TYPE_MATLAB '(rho0);']);
    
    % save the density matrix
    writeMatrix([pathname input_filename], rho0, 'rho0');

    % create grid (variables are used for interpolation)
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    
    % interpolate onto staggered grid in x direction and save
    rho0_sg = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x + kgrid.dx/2, kgrid.y, kgrid.z, '*linear');
    rho0_sg(isnan(rho0_sg)) = rho0(isnan(rho0_sg));
    writeMatrix([pathname input_filename], rho0_sg, 'rho0_sgx');
    
    % interpolate onto staggered grid in y direction and save
    rho0_sg = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y + kgrid.dy/2, kgrid.z, '*linear');
    rho0_sg(isnan(rho0_sg)) = rho0(isnan(rho0_sg));
    writeMatrix([pathname input_filename], rho0_sg, 'rho0_sgy');
    
    % interpolate onto staggered grid in z direction and save
    rho0_sg = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y, kgrid.z + kgrid.dz/2, '*linear');
    rho0_sg(isnan(rho0_sg)) = rho0(isnan(rho0_sg));
    writeMatrix([pathname input_filename], rho0_sg, 'rho0_sgz');
    
    % clear variable to free memory
    clear('kgrid', 'rho0', 'rho0_sg');
    
    % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

    % create the sound speed matrix
    c0       = c0_background * ones(Nx, Ny, Nz, 'single');   
    c0(ball) = c0_ball;
    
    % set the reference sound speed to the maximum in the medium
    c_ref = max(c0(:));

    % get the sound speed at the location of the source
    c_source = min(c0(:));
    
    % make sure the input is in the correct data format
    eval(['c0 = ' MATRIX_DATA_TYPE_MATLAB '(c0);']);
    
    % save the sound speed matrix
    writeMatrix([pathname input_filename], c0, 'c0');

    % clear variable to free memory
    clear('c0', 'ball');

    % :::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::---:::

    % make sure the inputs are in the correct data format
    eval(['alpha_coeff = ' MATRIX_DATA_TYPE_MATLAB '(alpha_coeff);']);
    eval(['alpha_power = ' MATRIX_DATA_TYPE_MATLAB '(alpha_power);']);

    % save the absorption variables
    writeMatrix([pathname input_filename], alpha_coeff, 'alpha_coeff');
    writeMatrix([pathname input_filename], alpha_power, 'alpha_power');

    % clear variables to free memory
    clear('alpha_coeff', 'alpha_power');

    % ---------------------------------------------------------------------
    % WRITE THE SOURCE PARAMETERS
    % ---------------------------------------------------------------------

    % update command line status
    toc; tic; fprintf('Writing source parameters... ');

    % define a square source mask facing in the x-direction using the
    % normal k-Wave syntax
    p_mask = false(Nx, Ny, Nz);
    p_mask(1 + pml_x_size, Ny/2 - source_y_size/2:Ny/2 + source_y_size/2, ...
        Nz/2 - source_z_size/2:Nz/2 + source_z_size/2) = 1;

    % find linear source indices
    p_source_index = find(p_mask == 1);
    p_source_index = reshape(p_source_index, [], 1);

    % make sure the input is in the correct data format
    eval(['p_source_index = ' INTEGER_DATA_TYPE_MATLAB '(p_source_index);']);

    % save the source index matrix
    writeMatrix([pathname input_filename], p_source_index, 'p_source_index');

    % clear variables to free memory
    clear p_mask p_source_index;

    % define a time varying sinusoidal source
    p_source_input  = source_strength .* ...
        sin(2 * pi * source_freq * (0:(Nt - 1)) * dt);

    % apply an cosine ramp to the beginning to avoid startup transients
    ramp_length = round((2 * pi / source_freq) / dt);
    p_source_input(1:ramp_length) = p_source_input(1:ramp_length) .* ...
        (-cos( (0:(ramp_length - 1)) * pi / ramp_length ) + 1)/2;

    % scale the source magnitude to be in the correct units for the code
    p_source_input = p_source_input .* (2 * dt ./ (3 * c_source * dx));
    
    % cast matrix to single precision
    eval(['p_source_input = ' MATRIX_DATA_TYPE_MATLAB '(p_source_input);']);

    % save the input signal
    writeMatrix([pathname input_filename], p_source_input, 'p_source_input');

    % clear variables to free memory
    clear('p_source_input');

    % ---------------------------------------------------------------------
    % WRITE THE SENSOR PARAMETERS
    % ---------------------------------------------------------------------

    % update command line status
    toc; tic; fprintf('Writing sensor parameters... ');

    % define a sensor mask through the central plane
    sensor_mask = false(Nx, Ny, Nz);
    sensor_mask(:, :, Nz/2) = 1;

    % extract the indices of the active sensor mask elements
    sensor_mask_index = find(sensor_mask);
    sensor_mask_index = reshape(sensor_mask_index, [], 1);

    % make sure the input is in the correct data format
    eval(['sensor_mask_index = ' INTEGER_DATA_TYPE_MATLAB '(sensor_mask_index);']);

    % save the sensor mask
    writeMatrix([pathname input_filename], sensor_mask_index, 'sensor_mask_index');

    % clear variables to free memory
    clear('sensor_mask', 'sensor_mask_index');

    % ---------------------------------------------------------------------
    % WRITE THE GRID PARAMETERS AND FILE ATTRIBUTES
    % ---------------------------------------------------------------------

    % update command line status
    toc; tic; fprintf('Writing grid parameters and attributes... ');

    % write grid parameters
    writeGrid([pathname input_filename], [Nx, Ny, Nz], [dx, dy, dz], ...
        [pml_x_size, pml_y_size, pml_z_size], ...
        [pml_x_alpha, pml_y_alpha, pml_z_alpha], ...
        Nt, dt, c_ref);

    % write flags
    writeFlags([pathname input_filename]);

    % set additional file attributes
    writeAttributes([pathname input_filename]);

    toc;
    
    % display the required syntax to run the C++ simulation
    disp(['Using a terminal window, navigate to the ' filesep 'binaries folder of the k-Wave Toolbox']);
    disp('Then, use the syntax shown below to run the simulation:');
    if isunix
        disp(['./kspaceFirstOrder-OMP -i ' pathname input_filename ' -o ' pathname output_filename ' --p_final --p_max']);
    else
        disp(['kspaceFirstOrder-OMP.exe -i "' pathname input_filename '" -o "' pathname output_filename '" --p_final --p_max']);
    end    

    return

% =========================================================================
% READ THE OUTPUT FILE AND PLOT VISUALISATION
% =========================================================================

else
    
    % load output data from the C++ simulation
    sensor_data.p_final = h5read([pathname output_filename], '/p_final');
    sensor_data.p_max   = h5read([pathname output_filename], '/p_max');

    % take an x-y slice through the final pressure output (this is recorded
    % over the entire grid)
    sensor_data.p_final = squeeze(sensor_data.p_final(:, :, Nz/2));

    % reshape the maximum pressure output (this is recorded at the grid points
    % specified by the sensor mask)
    sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Ny]);

    % add a display mask
    ball_outline = makeCircle(Nx, Ny, ball_x, ball_y, ball_radius);
    sensor_data.p_max  (ball_outline == 1) = max(sensor_data.p_max(:));
    sensor_data.p_final(ball_outline == 1) = max(sensor_data.p_final(:));

    % remove the pml
    sensor_data.p_max   = sensor_data.p_max  (1 + pml_x_size:end - pml_x_size, 1 + pml_y_size:end - pml_y_size);
    sensor_data.p_final = sensor_data.p_final(1 + pml_x_size:end - pml_x_size, 1 + pml_y_size:end - pml_y_size);

    % get a suitable plot scale
    x_vec = (0:(Nx - 2 * pml_x_size - 1)) * dx;
    y_vec = (0:(Ny - 2 * pml_y_size - 1)) * dy;
    [x_sc, scale, prefix] = scaleSI(max([x_vec, y_vec]));

    % plot the final pressure field in the x-y plane
    figure;
    subplot(1, 2, 1);
    imagesc(y_vec * scale, x_vec * scale, sensor_data.p_final, [-1, 1] * source_strength);
    colormap(getColorMap);
    xlabel(['y [' prefix 'm]']);
    ylabel(['x [' prefix 'm]']);
    axis image;
    title('Final Pressure Field');

    % plot the maximum pressure field in the x-y plane
    subplot(1, 2, 2);
    imagesc(y_vec * scale, x_vec * scale, sensor_data.p_max, [-2, 2] * source_strength);
    colormap(getColorMap);
    xlabel(['y [' prefix 'm]']);
    ylabel(['x [' prefix 'm]']);
    axis image;
    title('Maximum Pressure');

end