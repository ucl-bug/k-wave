% Focussed 2D Array With Directional Elements
%
% This example demonstrates the use of k-Wave to compute the outputs from a
% curved detector array which consists of several elements, each of which
% consists of a number of grid points.
%
% author: Ben Cox and Bradley Treeby
% date: 29th October 2010
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Ben Cox and Bradley Treeby

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
% SETUP THE GRID
% =========================================================================

% create the computational grid
Nx = 180;           % number of grid points in the x (row) direction
Ny = 180;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;    % [m/s]

% define the array of time points [s]
t_end = 12e-6;      % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE A FOCUSSED ARRAY OF DIRECTIONAL ELEMENTS
% =========================================================================

% define a semicircular sensor centred on the grid
semicircle_radius = 65; % [grid points]
arc = makeCircle(Nx, Ny, Nx/2, Ny/2, semicircle_radius, pi);

% find total number and indices of the grid points constituting the
% semicircle 
arc_indices = find(arc == 1);
Nv = length(arc_indices);

% calculate angles between grid points in the arc and the centre of the
% grid 
arc_angles = atan((kgrid.y(arc_indices)) ./ kgrid.x(arc_indices));

% sort the angles into ascending order, and adjust the indices accordingly
[sorted_arc_angles, sorted_index] = sort(arc_angles);
sorted_arc_indices = arc_indices(sorted_index);

% divide the semicircle into Ne separate sensor elements
Ne = 13;
sensor.mask = zeros(Nx, Ny);
for loop = 1:Ne
    
    % get the indices of the grid points belonging to the current element
    % (there is a two grid point gap between the elements)
    voxel_indices = sorted_arc_indices(floor((loop - 1) * Nv / Ne) + ...
        2:floor(loop * Nv / Ne) - 1);
    
    % add the element to the sensor.mask
    sensor.mask(voxel_indices) = 1;
    
end

% =========================================================================
% SIMULATION
% =========================================================================

% define an infinitely wide plane wave source (this is achieved by turning
% off the PML)
source.p_mask = zeros(Nx, Ny);
source.p_mask(140, :) = 1;

% set the display mask for the simulation
display_mask = source.p_mask + sensor.mask;

% define a time varying sinusoidal source
source_freq = 1e6;  % [Hz]
source_mag = 0.5;   % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
source.p = filterTimeSeries(kgrid, medium, source.p);

% run the simulation with the PML switched off on the sides, so that the
% source continues up to the edge of the domain (and from there infinitely,
% because of the periodic assumption implicit in pseudospectral methods)
input_args = {'PMLAlpha', [2, 0], 'DisplayMask', display_mask, 'PlotScale', [-0.75, 0.75]};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% split up the recorded data between the elements
element_data = zeros(Ne, kgrid.Nt);
for loop = 1:Ne
    
    % find indices of the sensor grid points in the sensor mask
    sensor_indices = find(sensor.mask == 1);
    
    % extract indices of the grid points belonging to the current element
    voxel_indices = sorted_arc_indices(floor((loop - 1) * Nv / Ne) + ...
        2:floor(loop * Nv / Ne) - 1);
    
    % indices of sensor_data that refer to the data for this element
    data_indices = zeros(length(voxel_indices), 1);
    for loop2 = 1:length(voxel_indices)
        data_indices(loop2) = find(sensor_indices == voxel_indices(loop2));
    end    

    % for one element per loop, average the time series from each of the
    % element's grid points to give one time series for the whole element
    element_data(loop, :) = mean(sensor_data(data_indices, :), 1);
    
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the elements of the array, and the source mask
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, display_mask);
colormap(flipud(gray));
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the time series recorded at each array element
figure;
stackedPlot(kgrid.t_array * 1e6, element_data);
xlabel('Time [\mus]');
ylabel('Time Series Recorded At Each Element');