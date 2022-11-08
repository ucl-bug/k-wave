function test_pass = pstdElastic2D_compare_labelled_and_binary_source_mask(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the simulation results using a labelled and
%     binary source mask
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 25th November 2014
%     last update - 19th February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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

% set additional literals to give further permutations of the test
COMPARISON_THRESH = 1e-15;
PML_INSIDE        = false;    
DATA_CAST         = 'off';

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 216;           % number of grid points in the x (row) direction
Ny = 216;           % number of grid points in the y (column) direction
dx = 50e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed_compression = 1500;  % [m/s]
medium.sound_speed_shear = 1000;        % [m/s]
medium.density = 1000;                  % [kg/m^3]

% create the time array
t_end = 20e-6;
kgrid.t_array = makeTime(kgrid, medium.sound_speed_compression, [], t_end);

% define a curved transducer element
arc_pos = [30, 30; 150, 30; 150, 200];
radius = [20, 30, 40];
diameter = [21, 15, 31];
focus_pos = [Nx/2, Ny/2];
[binary_mask, labelled_mask] = makeMultiArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);

% define a time varying sinusoidal source
source_freq = 0.25e6;       % [Hz]
source_mag = 0.5;           % [Pa]
source_1 = source_mag*sin(2*pi*source_freq*kgrid.t_array);

source_freq = 1e6;       % [Hz]
source_mag = 0.8;           % [Pa]
source_2 = source_mag*sin(2*pi*source_freq*kgrid.t_array);

source_freq = 0.05e6;       % [Hz]
source_mag = 0.2;           % [Pa]
source_3 = source_mag*sin(2*pi*source_freq*kgrid.t_array);

% assemble sources
labelled_sources(1, :) = source_1;
labelled_sources(2, :) = source_2;
labelled_sources(3, :) = source_3;

% assign sources for labelled source mask
source.s_mask = labelled_mask;
source.sxx = labelled_sources;
source.syy = labelled_sources;
source.sxy = labelled_sources;

% create a display mask to display the transducer
display_mask = binary_mask;

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1, 1, Nx, Ny].';

% set the record mode capture the final wave-field and the statistics at
% each sensor point 
sensor.record = {'p_final', 'p_max'};

% assign the input options
input_args = {'DataCast', DATA_CAST, 'DisplayMask', display_mask,...
    'PMLInside', PML_INSIDE, 'PlotPML', false, 'PlotSim', plot_simulations};

% run the simulation using the labelled source mask
sensor_data_labelled = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% reassign the source using a binary source mask
source.s_mask = binary_mask;
index_mask = labelled_mask(labelled_mask ~= 0);
source.sxx = labelled_sources(index_mask, :);
source.syy = source.sxx;
source.sxy = source.sxx;

% run the simulation using the a binary source mask
sensor_data_binary = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% compute the error from the first cuboid
L_inf_final = max(abs(sensor_data_labelled.p_final(:) - sensor_data_binary.p_final(:))) / max(abs(sensor_data_binary.p_final(:)));
L_inf_max   = max(abs(sensor_data_labelled.p_max(:)   - sensor_data_binary.p_max(:)))   / max(abs(sensor_data_binary.p_max(:)));

% compute pass
if (L_inf_max > COMPARISON_THRESH) || (L_inf_final > COMPARISON_THRESH)
    test_pass = false;
end

% ----------------------------------------

% repeat for velocity source
clear source;
source.u_mask = labelled_mask;
source.ux = labelled_sources * 1e-6;
source.uy = source.ux;

% run the simulation using the labelled source mask
sensor_data_labelled = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% reassign the source using a binary source mask
clear source;
source.u_mask = binary_mask;
index_mask = labelled_mask(labelled_mask ~= 0);
source.ux = labelled_sources(index_mask, :) * 1e-6;
source.uy = source.ux;

% run the simulation using the a binary source mask
sensor_data_binary = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% compute the error from the first cuboid
L_inf_final = max(abs(sensor_data_labelled.p_final(:) - sensor_data_binary.p_final(:))) / max(abs(sensor_data_binary.p_final(:)));
L_inf_max   = max(abs(sensor_data_labelled.p_max(:)   - sensor_data_binary.p_max(:)))   / max(abs(sensor_data_binary.p_max(:)));

% compute pass
if (L_inf_max > COMPARISON_THRESH) || (L_inf_final > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% PLOT COMPARISONS
% =========================================================================

if plot_comparisons

    % plot the final wave-field
    figure;
    subplot(1, 3, 1), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_labelled.p_final, [-1 1]);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('Binary');
    
    subplot(1, 3, 2), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, sensor_data_binary.p_final, [-1 1]);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('Labelled');

    subplot(1, 3, 3), imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, abs(sensor_data_binary.p_final - sensor_data_labelled.p_final), [-1 1]);
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
    title('Difference');

end