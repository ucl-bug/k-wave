% Defining A Source Using An Array Transducer Example
% 
% This example provides a demonstration of using the kWaveArray class to
% define an array transducer with three arc-shaped elements without
% staircasing errors.
%
% For a more detailed discussion of this example and the underlying
% techniques, see E. S. Wise, B. T. Cox, J. Jaros, & B. E. Treeby (2019).
% Representing arbitrary acoustic source and sensor distributions in
% Fourier collocation methods. The Journal of the Acoustical Society of
% America, 146(1), 278-288. https://doi.org/10.1121/1.5116132.
%
% author: Bradley Treeby
% date: 4th September 2018
% last update: 2nd November 2022
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2022 Bradley Treeby

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
% DEFINE KWAVEARRAY
% =========================================================================

% create empty array
karray = kWaveArray;

% define arc properties
radius    = 50e-3;              % [m]
diameter  = 30e-3;              % [m]
focus_pos = [-20e-3, 0];        % [m]

% add arc shaped element
elem_pos  = [10e-3, -40e-3];    % [m]
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [20e-3, 0];         % [m]
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [10e-3, 40e-3];     % [m]
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% move the array down 10 mm, and rotate by 10 degrees (this moves all the
% elements together)
karray.setArrayPosition([10e-3, 0], 10);

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 256;
dx = 0.5e-3;
Ny = 256;
dy = 0.5e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% medium properties
medium.sound_speed = 1500;

% time array
kgrid.makeTime(medium.sound_speed);

% =========================================================================
% SIMULATION
% =========================================================================

% assign binary mask from karray to the source mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% set source signals, one for each physical array element
f1 = 100e3;
f2 = 200e3;
f3 = 500e3;
sig1 = toneBurst(1/kgrid.dt, f1, 3);
sig2 = toneBurst(1/kgrid.dt, f2, 5);
sig3 = toneBurst(1/kgrid.dt, f3, 5);

% combine source signals into one array
source_signal = zeros(2, max(length(sig1), length(sig2)));
source_signal(1, 1:length(sig1)) = sig1;
source_signal(2, 1:length(sig2)) = sig2;
source_signal(3, 1:length(sig3)) = sig3;

% get distributed source signals (this automatically returns a weighted
% source signal for each grid point that forms part of the source)
source.p = karray.getDistributedSourceSignal(kgrid, source_signal);

% run k-Wave simulation (no sensor is used for this example)
kspaceFirstOrder2D(kgrid, medium, source, []);

% =========================================================================
% VISUALISATION
% =========================================================================

% create pml mask (default size in 2D is 20 grid points)
pml_size = 20;
pml_mask = false(Nx, Ny);
pml_mask(1:pml_size, :) = 1;
pml_mask(:, 1:pml_size) = 1;
pml_mask(end - pml_size + 1:end, :) = 1;
pml_mask(:, end - pml_size + 1:end) = 1;

% plot source and pml masks
figure;
imagesc(kgrid.y_vec, kgrid.x_vec, source.p_mask | pml_mask);
axis image;
colormap(flipud(gray));

% overlay the physical source positions
hold on;
karray.plotArray(false);
