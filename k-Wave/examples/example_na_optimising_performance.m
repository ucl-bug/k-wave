% Optimising k-Wave Performance Example
%
% This example demonstrates how to increase the computational performance
% of k-Wave using optional input parameters and data casting. A separate
% standardised benchmarking script benchmark is also included within the
% k-Wave toolbox to allow computational times to be compared across
% different computers and GPUs. 
%
% author: Bradley Treeby
% date: 15th July 2009
% last update: 6th August 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

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

example_number = 1;
% 1: default input options
% 2: nearest neighbour Cartesian interpolation and plotting switched off
% 3: as above with 'DataCast' set to 'single'
% 4: as aobve with 'DataCast' set to 'gpuArray-single'

% change scale to 2 to increase the computational time
scale = 1;

% =========================================================================
% SIMULATION
% =========================================================================

% assign the grid size and create the computational grid
Nx = 256 * scale;           % number of grid points in the x direction
Ny = 256 * scale;           % number of grid points in the y direction
x = 10e-3;                  % grid size in the x direction [m]
y = 10e-3;                  % grid size in the y direction [m]
dx = x / Nx;                % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% load the initial pressure distribution from an image and scale
p0_magnitude = 2;           % [Pa]
source.p0 = p0_magnitude * loadImage('EXAMPLE_source_two.bmp');
source.p0 = resize(source.p0, [Nx, Ny]);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% define a centered Cartesian circular sensor
sensor_radius = 4.5e-3;     % [m]
num_sensor_points = 100;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input options
switch example_number
    case 1
        
        % set input arguments
        input_args = {};
        
    case 2
        
        % set input arguments
        input_args = {'PlotSim', false};
        
        % convert Cartesian sensor mask to binary mask
        sensor.mask = cart2grid(kgrid, sensor.mask);
        
    case 3
        
        % set input arguments
        input_args = {'PlotSim', false, 'DataCast', 'single'};
        
        % convert Cartesian sensor mask to binary mask
        sensor.mask = cart2grid(kgrid, sensor.mask);
        
    case 4
        
        % set input arguments
        input_args = {'PlotSim', false, 'DataCast', 'gpuArray-single'};
        
        % convert Cartesian sensor mask to binary mask
        sensor.mask = cart2grid(kgrid, sensor.mask);
        
        % run GPU warmup
        tic; fft2(gpuArray(rand(128, 128, 'single'))); toc;
        tic; fft2(gpuArray(rand(128, 128, 'single'))); toc;
        tic; fft2(gpuArray(rand(128, 128, 'single'))); toc;
        
end

% run the simulation
kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});