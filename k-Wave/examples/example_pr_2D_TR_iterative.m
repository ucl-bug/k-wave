% 2D Iterative Image Improvement Using Time Reversal Example
%
% This example demonstrates how photoacoustic image reconstruction may be
% improved iteratively using time reversal and a positivity condition. This
% example builds on the 2D Time Reversal Reconstruction For A Line Sensor
% Example.
%
% author: Ben Cox, Bradley Treeby and Felix Lucka
% date: 22nd January 2012
% last update: 22nd October 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012- Ben Cox and Bradley Treeby

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
% SET UP AND RUN THE SIMULATION
% =========================================================================

% define literals
NUMBER_OF_ITERATIONS = 5;   % number of iterations
PML_SIZE = 20;              % size of the perfectly matched layer in grid points

% create the computational grid
Nx = 128 - 2 * PML_SIZE;    % number of grid points in the x direction
Ny = 256 - 2 * PML_SIZE;    % number of grid points in the y direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% create the time array
kgrid.makeTime(medium.sound_speed);

% load an image for the initial pressure distribution
p0_image = loadImage('EXAMPLE_k-Wave.png');

% make it binary
p0_image = double(p0_image > 0);

% smooth and scale the initial pressure distribution
p0 = smooth(p0_image, true);

% assign to the source structure
source.p0 = p0;

% define an L-shaped sensor mask
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;
sensor.mask(:, 1) = 1;

% set the input arguments: force the PML to be outside the computational
% grid, switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_SIZE, 'Smooth', false, ...
    'PlotPML', false, 'PlotSim', false}; 

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% RECONSTRUCT AN IMAGE USING ITERATIVE TIME REVERSAL
% =========================================================================

% set the initial reconstructed image to be zeros
p0_estimate = zeros(Nx, Ny);

% set the initial model times series to be zero
modelled_time_series = zeros(size(sensor_data));

% calculate the difference between the measured and modelled data
data_difference = sensor_data - modelled_time_series;

% track goodness of fit and relative error during the iteration
gof     = ones(NUMBER_OF_ITERATIONS+1, 1);
rel_err = ones(NUMBER_OF_ITERATIONS+1, 1); 

% =========================================================================
% ITERATE TO IMPROVE THE IMAGE
% =========================================================================

for loop = 1:NUMBER_OF_ITERATIONS
    
    % assign the data_difference as a time-reversal source
    source.p_mask = sensor.mask;
    source.p = fliplr(data_difference);
    source = rmfield(source, 'p0');
    source.p_mode = 'dirichlet';   

    % set the simulation to record the final image (at t = 0)
    sensor.record = {'p_final'};
    
    % run the time reversal reconstruction
    p0_update = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    % update the image    
    p0_estimate = p0_estimate + p0_update.p_final;
    
    % apply a positivity condition
    p0_estimate = p0_estimate .* (p0_estimate >= 0);
    
    % store the latest image estimate
    p0_iterates{loop} = p0_estimate;
    
    
    % remove the source used in the previous time reversal
    source = rmfield(source, 'p');

    % set the initial pressure to be the latest estimate of p0
    source.p0 = p0_estimate;
    
    % set the simulation to record the time series
    sensor = rmfield(sensor, 'record');
    
    % calculate the time series using the latest estimate of p0
    modelled_time_series = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    % calculate the error in the estimated time series
    data_difference = sensor_data - modelled_time_series;

    % measure goodness of fit and relative error
    gof(loop+1)     = norm(data_difference(:))^2/norm(sensor_data(:))^2;
    rel_err(loop+1) = norm(p0_estimate(:) - p0(:))^2/norm(p0(:))^2;
    
end

% =========================================================================
% VISUALISATION
% =========================================================================

% set the color scale
c_axis = [0, 1.1];

% plot the initial pressure
figure;
imagesc(p0, c_axis);
axis image;
set(gca, 'XTick', [], 'YTick', []);
title('Initial Acoustic Pressure');
colorbar;
scaleFig(1, 0.65);

% add lines indicating the sensor mask
hold on;
plot([Ny, 1], [1, 1], 'k-');
plot([1, 1], [1, Nx], 'k-');

% plot the first iteration
figure;
imagesc(p0_iterates{1}, c_axis);
axis image;
set(gca, 'XTick', [], 'YTick', []);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Time Reversal Reconstruction');
colorbar;
scaleFig(1, 0.65);

% add lines indicating the sensor mask
hold on;
plot([Ny, 1], [1, 1], 'k-');
plot([1, 1], [1, Nx], 'k-');

% plot the 2nd iteration
figure;
imagesc(p0_iterates{2}, c_axis);
axis image;
set(gca, 'XTick', [], 'YTick', []);
title('Time Reversal Reconstruction, 2 Iterations');
colorbar;
scaleFig(1, 0.65);

% add lines indicating the sensor mask
hold on;
plot([Ny, 1], [1, 1], 'k-');
plot([1, 1], [1, Nx], 'k-');

% plot the last iteration
figure;
imagesc(p0_iterates{end}, c_axis);
axis image;
set(gca, 'XTick', [], 'YTick', []);
title(['Time Reversal Reconstruction, ' int2str(NUMBER_OF_ITERATIONS) ' Iterations']);
colorbar;
scaleFig(1, 0.65);

% add lines indicating the sensor mask and 'visible region'
hold on;
plot([Ny, 1], [1, 1], 'k-');
plot([1, 1], [1, Nx], 'k-');
plot([Ny, 1], [1, Nx], 'k--');

% plot gof and relative error
figure();
plot(0:NUMBER_OF_ITERATIONS, gof, 'DisplayName', 'goodness of fit'); hold on
plot(0:NUMBER_OF_ITERATIONS, rel_err, 'DisplayName', 'relative error'); hold on
legend(gca,'show');
