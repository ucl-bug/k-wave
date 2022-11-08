% Defining A Sensor Using An Array Transducer Example
%
% This example provides a demonstration of using the kWaveArray class to
% define an array transducer with 20 arc-shaped elements which is then used
% as a receiver array. It builds on the Defining A Source Using An Array
% Transducer Example.
%
% author: Bradley Treeby
% date: 4th September 2018
% last update: 3rd November 2022
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

% set the properties for the arc shaped elements and the ring geometry in
% which they're placed
radius_of_curv = 100e-3;
diameter       = 8e-3;
ring_radius    = 50e-3;
num_elements   = 20;

% orient all elements towards the centre of the grid
focus_pos = [0, 0];

% generate the centre position for each element in Cartesian space using
% makeCartCircle (these positions could also be defined manually, etc)
elem_pos = makeCartCircle(ring_radius, num_elements, [0, 0]);

% add elements to the array
for ind = 1:num_elements
    karray.addArcElement(elem_pos(:, ind), radius_of_curv, diameter, focus_pos);
end

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
% CONVENTIONAL SIMULATION
% =========================================================================

% set source as a square (directional) and a disc
source.p0 = makeDisc(Nx, Ny, Nx/4 + 20, Ny/4, 3);
source.p0(100:120, 50:200) = 1;

% assign Cartesian points
sensor.mask = elem_pos;

% run the k-Wave simulation using point detectors in the normal way
sensor_data_point = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% KWAVEARRAY SIMULATION
% =========================================================================

% assign binary mask from karray to the sensor mask
sensor.mask = karray.getArrayBinaryMask(kgrid);

% run k-Wave simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% combine data to give one trace per physical array element
combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);

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

% plot source, sensor, and pml masks
figure;
imagesc(kgrid.y_vec, kgrid.x_vec, sensor.mask | source.p0 | pml_mask);
axis image;
colormap(flipud(gray));

% overlay the physical source positions
hold on;
karray.plotArray(false);

% plot recorded sensor data
figure;
subplot(2, 1, 1)
imagesc(kgrid.t_array * 1e6, 1:num_elements, sensor_data_point);
xlabel('Time [\mus]');
ylabel('Detector Number');
title('Cartesian point detectors');
colorbar;

subplot(2, 1, 2);
imagesc(kgrid.t_array * 1e6, 1:num_elements, combined_sensor_data);
xlabel('Time [\mus]');
ylabel('Detector Number');
title('Arc detectors');
colorbar;

% plot a trace from the recorded sensor data
figure;
plot(kgrid.t_array * 1e6, sensor_data_point(1, :));
hold on;
plot(kgrid.t_array * 1e6, combined_sensor_data(1, :));
xlabel('Time [\mus]');
ylabel('Pressure [pa]');
legend('Cartesian point detectors', 'Arc detectors', 'Location', 'best');