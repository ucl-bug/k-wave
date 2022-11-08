function test_pass = kspaceFirstOrder3D_stream_to_disk(~, plot_simulations)
% DESCRIPTION:
%     Unit test to compare the output of kspaceFirstOrder3D with and
%     without the option to stream to disk.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 18th August 2014
%     last update - 7th August 2017
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
    plot_simulations = true;
end

% set pass variable
test_pass = true;

% set additional literals to give further permutations of the test
COMPARISON_THRESH = 1e-15;
DATA_CAST         = 'off';

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
dz = 0.1e-3;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny, Nz);	% [m/s]
medium.sound_speed(1:Nx/2, :, :) = 1800;    % [m/s]
medium.density = 1000*ones(Nx, Ny, Nz);     % [kg/m^3]
medium.density(:, Ny/4:end, :) = 1200;      % [kg/m^3]

% create initial pressure distribution using makeBall
ball_magnitude = 10;    % [Pa]
ball_x_pos = 38;        % [grid points]
ball_y_pos = 32;        % [grid points]
ball_z_pos = 32;        % [grid points]
ball_radius = 5;        % [grid points]
ball_1 = ball_magnitude*makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

ball_magnitude = 10;    % [Pa]
ball_x_pos = 20;        % [grid points]
ball_y_pos = 20;        % [grid points]
ball_z_pos = 20;        % [grid points]
ball_radius = 3;        % [grid points]
ball_2 = ball_magnitude*makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

source.p0 = ball_1 + ball_2;

% define a binary sensor mask
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, :, Nz/4) = 1;

% input arguments
input_args = {'PlotSim', plot_simulations, 'DataCast', DATA_CAST};

% run the simulations
sensor_data_ref      = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_data_stream_1 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'StreamToDisk', true);
sensor_data_stream_2 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'StreamToDisk', 30);

% compute errors
L_inf_1 = max(abs(sensor_data_stream_1(:) - sensor_data_ref(:))) / max(abs(sensor_data_ref(:)));
L_inf_2 = max(abs(sensor_data_stream_2(:) - sensor_data_ref(:))) / max(abs(sensor_data_ref(:)));

% add a start index
sensor.record_start_index = 50;

% run the simulations
sensor_data_ref      = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sensor_data_stream_1 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'StreamToDisk', true);
sensor_data_stream_2 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'StreamToDisk', 30);

% compute errors
L_inf_3 = max(abs(sensor_data_stream_1(:) - sensor_data_ref(:))) / max(abs(sensor_data_ref(:)));
L_inf_4 = max(abs(sensor_data_stream_2(:) - sensor_data_ref(:))) / max(abs(sensor_data_ref(:)));

% get maximum error
L_inf_max = max([L_inf_1, L_inf_2, L_inf_3, L_inf_4]);

% compute pass
if (L_inf_max > COMPARISON_THRESH)
    test_pass = false;
end