% Snell's Law And Critical Angle Reflection Example
%
% This example illustrates Snell's law by steering a tone burst from a
% linear array transducer within a layered heterogeneous medium. It builds
% on the Simulating Transducer Field Patterns Example.
%
% author: Bradley Treeby
% date: 10th December 2009
% last update: 4th May 2017
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

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = Nx;            % number of grid points in the y (column) direction
dx = 50e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define sound speed in the two layers
c1 = 1500;          % [m/s]
c2 = 3000;          % [m/s]

% define the properties of the propagation medium
medium.sound_speed = c1 * ones(Nx, Ny);     % [m/s]
medium.sound_speed(Nx/2:end, :) = c2;       % [m/s]
medium.density = 1000;                      % [kg/m^3]
medium.alpha_coeff = 0.75;                  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% create the time array
kgrid.makeTime(medium.sound_speed);

% define a source mask for a linear element transducer with an odd number
% of elements 
num_elements = 61;
source.p_mask = zeros(Nx, Ny);
x_offset = 25;
y_offset = 20;
start_index = Ny/2 - round(num_elements/2) + 1 - y_offset;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% define the properties of the tone burst
sampling_freq = 1/kgrid.dt;     % [Hz]
steering_angle = 35;            % [deg]
element_spacing = dx;           % [m]
tone_burst_freq = 1e6;          % [Hz]
tone_burst_cycles = 8;

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

% use geometric beam forming to calculate the tone burst offsets for each
% transducer element based on the element index
tone_burst_offset = 200 + element_spacing * element_index * ...
    sin(steering_angle * pi/180) / (c1 * kgrid.dt);

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
    'SignalOffset', tone_burst_offset);

% assign the input options
input_args = {'DisplayMask', source.p_mask};

% run the simulation
kspaceFirstOrder2D(kgrid, medium, source, [], input_args{:});