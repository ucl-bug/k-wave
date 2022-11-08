% Automatic Sound Speed Selection Example
%
% This example illustrates an autofocus approach for automatically
% selecting the sound speed used in photoacoustic image reconstruction. The
% selection is based on maximising the sharpness of the reconstructed image
% as quantified by a focus function.
%
% For a more detailed discussion of this example and the underlying
% techniques, see B. E. Treeby, T. K. Varslot, E. Z. Zhang, J. G. Laufer,
% and P. C. Beard, "Automatic sound speed selection in photoacoustic image
% reconstruction using an autofocus approach," J. Biomed. Optics, vol. 16,
% no. 9, p. 090501, 2011.    
%
% author: Bradley Treeby
% date: 26th January 2012
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2019 Bradley Treeby

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

% load the initial pressure distribution from an image and scale
p0_magnitude = 2;
p0 = p0_magnitude * loadImage('EXAMPLE_source_two.bmp');

% assign the grid size and create the computational grid
PML_size = 20;              % size of the PML in grid points
Nx = 256 - 2 * PML_size;    % number of grid points in the x direction
Ny = 256 - 2 * PML_size;    % number of grid points in the y direction
x = 10e-3;                  % total grid size [m]
y = 10e-3;                  % total grid size [m]
dx = x / Nx;                % grid point spacing in the x direction [m]
dy = y / Ny;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% resize the input image to the desired number of grid points
p0 = resize(p0, [Nx, Ny]);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(p0, true);

% assign to the source structure
source.p0 = p0;

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% define a centered circular sensor mask
radius = 100;
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, radius);

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input options
input_args = {'Smooth', false, 'PMLInside', false, 'PlotPML', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% RECONSTRUCTION
% =========================================================================

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the reconstruction in a loop using different values for sound speed
sound_speed_array = 1450:5:1550;
focus_func = zeros(3, length(sound_speed_array));
for index = 1:length(sound_speed_array)
    
    % update value of sound speed
    medium.sound_speed = sound_speed_array(index);
    
    % run the reconstruction
    p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    % update the value of the focus function
    focus_func(1, index) = sharpness(p0_recon, 'Brenner');
    focus_func(2, index) = sharpness(p0_recon, 'Tenenbaum');
    focus_func(3, index) = sharpness(p0_recon, 'NormVariance');
    
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the variation in the focus function with sound speed
figure;
plot(sound_speed_array, focus_func(1, :) ./ max(focus_func(1, :)), 'k-');
hold on;
plot(sound_speed_array, focus_func(2, :) ./ max(focus_func(2, :)), 'r-');
plot(sound_speed_array, focus_func(3, :) ./ max(focus_func(3, :)), 'b-');
ylabel('Normalised Focus Function');
xlabel('Sound Speed [m/s]');
legend('Brenner Gradient', 'Tenenbaum Gradient', 'Normalised Variance');