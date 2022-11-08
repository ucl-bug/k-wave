% Saving Movie Files Example
%
% This example demonstrates how to save the simulation animations as a
% movie. It builds on the Heterogeneous Propagation Medium Example.
%
% author: Bradley Treeby
% date: 30th June 2009
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

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium    
medium.sound_speed = 1500 * ones(Nx, Ny); % [m/s]
medium.sound_speed(1:Nx/2, :) = 1800;   % [m/s]
medium.density = 1000 * ones(Nx, Ny);     % [kg/m^3]
medium.density(:, Ny/4:end) = 1200;     % [kg/m^3]

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2;

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% define the input arguments
input_args = {'RecordMovie', true, 'MovieName', 'example_movie_1'};

% run the simulation
kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% get allowable VideoWriter profiles
profiles = VideoWriter.getProfiles();

% if MPEG profile exists, run second simulation
if any(ismember({profiles.Name}, 'MPEG-4'))

    % define a second set of input arguments
    input_args = {'MeshPlot', true, 'DisplayMask', 'off', 'PlotPML', false, ...
        'PlotFreq', 5, 'RecordMovie', true, 'MovieName', 'example_movie_2', ...
        'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};

    % run the simulation again
    kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

else
    
    % display warning about MPEG profile
    disp('The optional input ''MovieProfile'', ''MPEG-4'' is not supported in your version of MATLAB');
    
end