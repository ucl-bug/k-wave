function test_pass = pstdElastic2D_save_movie(~, ~)
% DESCRIPTION:
%     Unit test to generate movies using pstdElastic2D.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 19th February 2017
%     last update - 18th June 2017
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

% set pass variable
test_pass = true;

% folder to store temporary movies in
if nargin == 0
    folder = '';
else
    folder = tempdir;
end

% movie names
movie_name_1 = 'test-movie-elastic-2D-avi';
movie_name_2 = 'test-movie-elastic-2D-mpg';

% =========================================================================
% DEFAULT SETTINGS
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500*ones(Nx, Ny); % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);     % [m/s]
medium.density                 = 1000*ones(Nx, Ny); % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000; % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 800;  % [m/s]
medium.density(Nx/2:end, :)                 = 1200; % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.5; % [dB/(MHz^2 cm)]

% create the time array
cfl   = 0.1;        % Courant-Friedrichs-Lewy number
t_end = 8e-6;       % [s]
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed_compression(:)), cfl, t_end);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 30;    % [grid points]
disc_y_pos = 64;    % [grid points]
disc_radius = 5;    % [grid points]
source.p0 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define a centered circular sensor
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, 20);

% define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
display_mask(Nx/2 - 1, :) = 1;

% define input arguments
input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
    'DisplayMask', display_mask, 'DataCast', 'single', ...
    'RecordMovie', true, 'MovieName', [folder movie_name_1]};

% run the simulation
pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% check movie file was created, if so, delete it
if exist([folder movie_name_1 '.avi'], 'file')
    if nargin ~= 0
        delete([folder movie_name_1 '.avi']);
    end
else
    test_pass = false;
end

% =========================================================================
% CUSTOM SETTINGS
% =========================================================================

% get list of profiles
profiles = VideoWriter.getProfiles;

% check if profile is supported, if so, run test
if any(ismember({profiles.Name}, 'MPEG-4'))

    % define input arguments
    input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
        'DisplayMask', display_mask, 'DataCast', 'single', ...
        'RecordMovie', true, 'MovieName', [folder movie_name_2], ...
        'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};

    % run the simulation
    pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

    % check movie file was created, if so, delete it
    if exist([folder movie_name_2 '.mp4'], 'file')
        if nargin ~= 0
            delete([folder movie_name_2 '.mp4']);
        end
    else
        test_pass = false;
    end
    
end