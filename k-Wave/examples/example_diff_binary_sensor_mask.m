% Using A Binary Sensor Mask Example
%
% This example provides a simple demonstration of using k-Wave to model the
% diffusion of heat within a two-dimensional homogeneous medium with a
% constant source term (volume rate of heat deposition). It uses a binary
% source mask to extract the temperature time profile. It builds on the
% Heat Diffusion In A Homogeneous Medium and Constant Rate Of Heat
% Deposition examples.
%
% author: Bradley Treeby
% date: 15th January 2017
% last update: 28th April 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby

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
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 1e-3;          % grid point spacing in the x direction [m]
dy = 1e-3;          % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define medium properties
medium.density              = 1079;     % [kg/m^3]
medium.thermal_conductivity = 0.52;     % [W/(m.K)]
medium.specific_heat        = 3540;     % [J/(kg.K)]

% set initial temperature distribution to be constant
source.T0 = 37;

% set Gaussian volume rate of heat deposition
width = 4*dx;
source.Q = 2e6 .* exp( -(kgrid.x ./ width).^2 - (kgrid.y ./ width).^2 );

% create binary sensor mask with a single sensor point in the centre of the
% grid
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, Ny/2) = 1;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, sensor);

% set source on time and off time
on_time  = 60;    % [s]
off_time = 120;   % [s]

% set time step size
dt = 0.1;

% take time steps
kdiff.takeTimeStep(round(on_time / dt), dt);

% turn off heat source
kdiff.Q = 0;

% take time steps
kdiff.takeTimeStep(round(off_time / dt), dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% create time axis
t_axis = (0:kdiff.time_steps_taken - 1) * dt;

% plot temperature time profile
figure;
plot(t_axis, kdiff.sensor_data);
xlabel('Time [s]');
ylabel('Temperature [^\circC]');