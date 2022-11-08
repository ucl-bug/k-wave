function test_pass = kspaceFirstOrder1D_save_movie(~, ~)
% DESCRIPTION:
%       Unit test to generate movies using kspaceFirstOrder1D.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 19th February 2017
%       last update - 19th February 2017
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
movie_name_1 = 'test-movie-1D-avi';
movie_name_2 = 'test-movie-1D-mpg';

% =========================================================================
% DEFAULT SETTINGS
% =========================================================================

% create the computational grid
Nx = 512;       % number of grid points in the x (row) direction
dx = 0.05e-3;   % grid point spacing in the x direction [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, 1);    	% [m/s]
medium.sound_speed(1:round(Nx/3)) = 2000;	% [m/s]
medium.density = 1000*ones(Nx, 1);          % [kg/m^3]
medium.density(round(4*Nx/5):end) = 1500;   % [kg/m^3]

% create initial pressure distribution using a smoothly shaped sinusoid
x_pos = 280;    % [grid points]
width = 100;    % [grid points]
height = 1;     % [au]
in = (0:pi/(width/2):2*pi).';
source.p0 = [zeros(x_pos, 1); ((height/2)*sin(in-pi/2)+(height/2)); zeros(Nx - x_pos  - width - 1, 1)];

% create a Cartesian sensor mask
sensor.mask = [-10e-3, 10e-3];  % [mm]

% define the input arguments
input_args = {'PlotPML', false, 'RecordMovie', true, 'MovieName', [folder movie_name_1]};

% run the simulation
kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

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

% get allowable VideoWriter profiles
profiles = VideoWriter.getProfiles();

% if MPEG profile exists, run second test
if any(ismember({profiles.Name}, 'MPEG-4'))

    % define the input arguments
    input_args = {'PlotPML', false, 'RecordMovie', true, 'MovieName', [folder movie_name_2], 'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};

    % run the simulation
    kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

    % check movie file was created, if so, delete it
    if exist([folder movie_name_2 '.mp4'], 'file')
        if nargin ~= 0
            delete([folder movie_name_2 '.mp4']);
        end
    else
        test_pass = false;
    end

end