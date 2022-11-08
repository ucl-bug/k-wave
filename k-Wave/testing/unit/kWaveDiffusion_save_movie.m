function test_pass = kWaveDiffusion_save_movie(~, ~)
% DESCRIPTION:
%     Unit test to generate movies using the kWaveDiffusion class.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 17th August 2015
%     last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Bradley Treeby

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

%#ok<*INUSD>

% set pass variable
test_pass = true;

% folder to store temporary movies in
folder = tempdir;

% movie names
movie_name_1D = 'test-movie-1D.avi';
movie_name_2D = 'test-movie-2D.avi';
movie_name_3D = 'test-movie-3D.avi';
movie_name_1D_mpg = 'test-movie-1D-MPEG.mp4';
movie_name_2D_mpg = 'test-movie-2D-MPEG.mp4';
movie_name_3D_mpg = 'test-movie-3D-MPEG.mp4';

% get allowable VideoWriter profiles
profiles = VideoWriter.getProfiles();

% check if MPG profile exists
run_mpg_tests = any(ismember({profiles.Name}, 'MPEG-4'));

% =========================================================================
% 1D SIMULATION
% =========================================================================

% create grid
Nx = 128;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% define medium properties
medium.density              = 1079;     % [kg/m^3]
medium.thermal_conductivity = 0.52;     % [W/(m.K)]
medium.specific_heat        = 3540;     % [J/(kg.K)]

% set Gaussian initial temperature distribution
offset = 0;
width = dx;
source.T0 = exp(-((kgrid.x-offset)/width).^2) + 37;

% set input args
input_args = {'PlotScale', [37, 38], 'PlotFreq', 5,...
    'RecordMovie', true};

% set input args for saving MPG
input_args_mpg = [input_args, {'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}}];

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:}, 'MovieName', [folder movie_name_1D]);

% take time steps
Nt = 300;
dt = 0.5;
kdiff.takeTimeStep(Nt, dt);

% if MPEG profile exists, run second test
if run_mpg_tests

    % create kWaveDiffusion object:
    kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args_mpg{:}, 'MovieName', [folder movie_name_1D_mpg]);

    % take time steps
    kdiff.takeTimeStep(Nt, dt);

end

% =========================================================================
% 2D SIMULATION
% =========================================================================

% create grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% set Gaussian initial temperature distribution
source.T0 = makeDisc(Nx, Ny, Nx/2, Ny/2, 10);
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:}, 'MovieName', [folder movie_name_2D]);

% take time steps
kdiff.takeTimeStep(Nt, dt);

% if MPEG profile exists, run second test
if run_mpg_tests

    % create kWaveDiffusion object:
    kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args_mpg{:}, 'MovieName', [folder movie_name_2D_mpg]);

    % take time steps
    kdiff.takeTimeStep(Nt, dt);

end

% =========================================================================
% 3D SIMULATION
% =========================================================================

% create grid
Nx = 64;
Ny = 64;
Nz = 64;
dx = 1e-3;
dy = 1e-3;
dz = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% set Gaussian initial temperature distribution
source.T0 = makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, 10);
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:}, 'MovieName', [folder movie_name_3D]);

% take time steps
kdiff.takeTimeStep(Nt, dt);

% if MPEG profile exists, run second test
if run_mpg_tests

    % create kWaveDiffusion object:
    kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args_mpg{:}, 'MovieName', [folder movie_name_3D_mpg]);

    % take time steps
    kdiff.takeTimeStep(Nt, dt);

end

% =========================================================================
% DELETE FILES
% =========================================================================

% check movie files were created, if so, delete them
if exist([folder movie_name_1D], 'file')
    if nargin ~= 0
        delete([folder movie_name_1D]);
    end
else
    test_pass = false;
end
if exist([folder movie_name_2D], 'file')
    if nargin ~= 0
        delete([folder movie_name_2D]);
    end
else
    test_pass = false;
end
if exist([folder movie_name_3D], 'file')
    if nargin ~= 0
        delete([folder movie_name_3D]);
    end
else
    test_pass = false;
end

% check mpg movie files were created, if so, delete them
if run_mpg_tests
    if exist([folder movie_name_2D_mpg], 'file')
        if nargin ~= 0
            delete([folder movie_name_2D_mpg]);
        end
    else
        test_pass = false;
    end
    if exist([folder movie_name_1D_mpg], 'file')
        if nargin ~= 0
            delete([folder movie_name_1D_mpg]);
        end
    else
        test_pass = false;
    end
    if exist([folder movie_name_3D_mpg], 'file')
        if nargin ~= 0
            delete([folder movie_name_3D_mpg]);
        end
    else
        test_pass = false;
    end
end