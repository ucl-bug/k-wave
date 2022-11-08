% Setting An Initial Pressure Gradient Example
%
% This example demonstrates how to set an initial pressure gradient using
% kspaceSecondOrder. It builds on the Comparison Of Modelling Functions
% Example.
%
% author: Bradley Treeby
% date: 28th October 2010
% last update: 31st July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby

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

%#ok<*UNRCH> 

clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x (row) direction
Ny = 64;            % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_power = 1.5;   % [dB/(MHz^y cm)]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]

% set plot_frames to true to produce the plots given in the documentation
% or false to visualise the propagation of the pressure field
plot_frames = true;
if plot_frames
    
    % calculate and plot the pressure at a particular value of t
    Nt = 2;
    dt = 1e-6;  % [s]
    kgrid.setTime(Nt, dt);
    
else
    
    % calculate the plot the progression of the pressure field with t
    dt = 5e-9;  % [s]
    Nt = round(5e-6 / dt);
    kgrid.setTime(Nt, dt);
    
end

% create source distribution using makeDisc
disc_magnitude = 4; % [Pa]
disc_x_pos = 25;    % [grid points]
disc_y_pos = 40;    % [grid points]
disc_radius = 4;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 40;    % [grid points]
disc_y_pos = 25;    % [grid points]
disc_radius = 3;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source_distribution = disc_1 + disc_2;

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, Ny/2) = 1;

% define the input arguments
input_args = {'PlotFrames', plot_frames, 'MeshPlot', true, ...
    'PlotScale', [0 3], 'ExpandGrid', true};

% assign the source distribution to the initial pressure
source.p0 = source_distribution;
kspaceSecondOrder(kgrid, medium, source, sensor, input_args{:});

% assign the source distribution to the initial pressure gradient
source = rmfield(source, 'p0');
source.dp0dt = 5e6 * source_distribution;
kspaceSecondOrder(kgrid, medium, source, sensor, input_args{:});

% assign the source distribution to both the initial pressure and the
% initial pressure gradient 
source.p0 = source_distribution;
kspaceSecondOrder(kgrid, medium, source, sensor, input_args{:});