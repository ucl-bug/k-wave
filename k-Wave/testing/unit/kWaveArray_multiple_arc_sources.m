function test_pass = kWaveArray_multiple_arc_sources(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to test the generation of sources using kWaveArray.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 6th September 2018
%       last update - 4th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set pass variable
test_pass = true;

% set comparison threshold
comparison_thresh = 1e-13;

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

% options for offgrid
offgrid_options = {...
    'BLITolerance', 0.1, ...
    'BLIType', 'sinc', ...
    'UpsamplingRate', 8};

% create empty array
karray = kWaveArray(offgrid_options{:});

% add arc shaped element
elem_pos  = [30e-3, 0];
radius    = 50e-3;
diameter  = 30e-3;
focus_pos = [0, 1];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [-30e-3, 0];
radius    = 50e-3;
diameter  = 30e-3;
focus_pos = [0, 0];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [-50e-3, -50e-3];
radius    = 80e-3;
diameter  = 20e-3;
focus_pos = [0, 0];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% get grid weights for each source for manual comparison
source_weights_1 = karray.getElementGridWeights(kgrid, 1);
source_weights_2 = karray.getElementGridWeights(kgrid, 2);
source_weights_3 = karray.getElementGridWeights(kgrid, 3);

% get mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% plot
if plot_comparisons
    figure;
    imagesc(source.p_mask);
    axis image;
    colormap(flipud(gray));
end

% set source signals
f1 = 100e3;
f2 = 200e3;
f3 = 500e3;
sig1 = toneBurst(1/kgrid.dt, f1, 3);
sig2 = toneBurst(1/kgrid.dt, f2, 5);
sig3 = toneBurst(1/kgrid.dt, f3, 5);
source_signal = zeros(2, max(length(sig1), length(sig2)));
source_signal(1, 1:length(sig1)) = sig1;
source_signal(2, 1:length(sig2)) = sig2;
source_signal(3, 1:length(sig3)) = sig3;

% get distributed source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_signal);

% set a circular sensor mask
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, Nx/2 - 2);

% set input arguments
input_args = {...
    'PlotSim', plot_simulations, ...
    'PMLSize', 'auto', ...
    'PMLInside', false};

% run k-Wave simulation with combined sources
sensor_data_karray = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% run individual simulations directly
clear source;
source.p_mask = (source_weights_1 ~= 0);
source.p = source_weights_1(source_weights_1 ~= 0) * sig1;
sensor_data_1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

clear source;
source.p_mask = (source_weights_2 ~= 0);
source.p = source_weights_2(source_weights_2 ~= 0) * sig2;
sensor_data_2 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

clear source;
source.p_mask = (source_weights_3 ~= 0);
source.p = source_weights_3(source_weights_3 ~= 0) * sig3;
sensor_data_3 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% sum individual simulations
sensor_data_ind = sensor_data_1 + sensor_data_2 + sensor_data_3;

% compute error
err = max(sensor_data_karray(:) - sensor_data_ind(:));

% check values are the same
if err > comparison_thresh
    test_pass = false;
end

if plot_comparisons
   
    figure;
    subplot(3, 1, 1);
    imagesc(sensor_data_karray);
    title('karray');
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(sensor_data_ind);
    title('individual');
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(abs(sensor_data_karray - sensor_data_ind));
    title('difference');
    colorbar;
    
end