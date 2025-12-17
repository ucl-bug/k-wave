function test_pass = hdf5_save_in_parts(~, ~)
% DESCRIPTION:
%     Unit test to test to check the grid parameters that are written to
%     disk when saving in parts with the parameters written by
%     kspaceFirstOrder3D. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th June 2017
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017- Bradley Treeby

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

% pathname for the input file
pathname = tempdir;

% input filename (this must have the .h5 extension)
filename_1  = 'example_input_1.h5';
filename_2  = 'example_input_2.h5';

% loop twice, once for even grid sizes, once for odd
for ind = 0:1

    % delete files in case they already exist
    if exist([pathname filename_1], 'file')
        delete([pathname, filename_1]);
    end
    if exist([pathname filename_2], 'file')
        delete([pathname, filename_2]);
    end

    % =====================================================================
    % K-WAVE SIMULATION
    % =====================================================================

    % create the computational grid
    Nx = 60 + ind; 
    Ny = 62 + ind;
    Nz = 64 + ind;
    dx = 0.1e-3;
    dy = 0.1e-3;
    dz = 0.1e-3;
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % set pml size
    pml_size = 8;
    pml_alpha = 1.5;

    % define the properties of the propagation medium
    medium.sound_speed_ref = 1400;
    medium.sound_speed = 1500 * ones(Nx, Ny, Nz);
    medium.density = 1000 * ones(Nx, Ny, Nz);

    % create time array
    kgrid.makeTime(medium.sound_speed);

    % create initial pressure distribution using makeBall
    source.p0 = 3.4 * makeBall(Nx, Ny, Nz, 38, 32, 32, 5);

    % define sensor mask
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(round(Nx/2), round(Ny/2), :) = 1;

    % set input options
    input_args = {'SaveToDisk', [pathname, filename_1], 'PMLSize', pml_size, ...
        'PMLAlpha', pml_alpha, 'Smooth', false};

    % save to disk
    kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

    % =====================================================================
    % SAVE IN PARTS
    % =====================================================================

    % save sound speed
    writeMatrix([pathname, filename_2], single(medium.sound_speed), 'c0');

    % save density
    writeMatrix([pathname, filename_2], single(medium.density), 'rho0');
    writeMatrix([pathname, filename_2], single(medium.density), 'rho0_sgx');
    writeMatrix([pathname, filename_2], single(medium.density), 'rho0_sgy');
    writeMatrix([pathname, filename_2], single(medium.density), 'rho0_sgz');

    % save source
    writeMatrix([pathname, filename_2], single(source.p0), 'p0_source_input');

    % save sensor
    sensor_mask_index = find(sensor.mask);
    sensor_mask_index = reshape(sensor_mask_index, [], 1);
    writeMatrix([pathname, filename_2], uint64(sensor_mask_index), 'sensor_mask_index');

    % save grid parameters
    writeGrid([pathname, filename_2], [Nx, Ny, Nz], [dx, dy, dz], ...
        [pml_size, pml_size, pml_size], ...
        [pml_alpha, pml_alpha, pml_alpha], ...
        kgrid.Nt, kgrid.dt, medium.sound_speed_ref);

    % write flags
    writeFlags([pathname, filename_2]);

    % set additional file attributes
    writeAttributes([pathname, filename_2]);

    % =====================================================================
    % COMPARE FILES
    % =====================================================================

    % check input files were created
    if ~exist([pathname filename_1], 'file')
        test_pass = false;
    end
    if ~exist([pathname filename_2], 'file')
        test_pass = false;
    end

    % check input files are the same (there will be one difference because
    % of the file date)
    if h5compare([pathname, filename_1], [pathname, filename_2]) > 1
        test_pass = false;
    end

    % delete files
    delete([pathname, filename_1]);
    delete([pathname, filename_2]);
    
end