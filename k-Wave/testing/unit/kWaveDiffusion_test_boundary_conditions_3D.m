function test_pass = kWaveDiffusion_test_boundary_conditions_3D(plot_comparisons, plot_simulations) 
% DESCRIPTION:
%       Unit test to compare the output from the k-Wave diffusion class
%       for insulating and conducting boundary conditions. Requires the DTT
%       functions.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 9th April 2019
%       last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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

% check for dtt functions
if ~exist('dtt1D', 'file') 
    disp('WARNING: kWaveDiffusion_test_boundary_conditions_2D not tested, as dtt library is not present on path');
    return
end

% set comparison thresholds
insulating_heterog_diff = 2e-4;
insulating_homog        = 1e-14;
conducting_max_temp     = 2e-2;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% =========================================================================
% MATERIAL PROPERTIES
% =========================================================================

% define medium properties
medium.density              = 1079;     % [kg/m^3]
medium.thermal_conductivity = 0.52;     % [W/(m.K)]
medium.specific_heat        = 3540;     % [J/(kg.K)]

% define medium properties with perfusion
medium_perfused = medium;

% define blood properties
blood_density             = 1060;	% [kg/m^3]
blood_specific_heat       = 3617;	% [J/(kg.K)]
blood_perfusion_rate      = 0.01;  	% [1/s]

% define perfusion coefficient
medium_perfused.perfusion_coeff = blood_density .* blood_perfusion_rate .* blood_specific_heat ./ (medium.density .* medium.specific_heat);

% define arterial temperature
medium_perfused.blood_ambient_temperature = 0;  % [degC]

% set input args
input_args = {'PlotSim', plot_simulations};

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% create grid
Nx = 65;
Ny = 65;
Nz = 65;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define initial temperature distribution
source.T0 = 25 * makeBall(Nx, Ny, Nz, 0, 0, 0, 5);
source.T0 = smooth(source.T0, true);

% define homogeneous temperature distribution
source_homog.T0 = ones(Nx, Ny, Nz);

% =========================================================================
% INSULATING - HOMOGENEOUS TEMP
% =========================================================================

% set boundary condition
medium.boundary_condition = 'insulating';

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source_homog, [], input_args{:});

% time steps
Nt = 20;
dt = 1;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
temp_diff = max(kdiff.T(:)) - min(kdiff.T(:));

% check for test pass
if (temp_diff > insulating_homog)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end

% =========================================================================
% INSULATING - HOMOGENEOUS TEMP, PERFUSION
% =========================================================================

% set boundary condition
medium_perfused.boundary_condition = 'insulating';

% define arterial temperature
medium_perfused.blood_ambient_temperature = 1;  % [degC]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium_perfused, source_homog, [], input_args{:});

% time steps
Nt = 20;
dt = 1;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
temp_diff = max(kdiff.T(:)) - min(kdiff.T(:));

% check for test pass
if (temp_diff > insulating_homog)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end

% =========================================================================
% INSULATING - HETEROGENEOUS TEMP
% =========================================================================

% set boundary condition
medium.boundary_condition = 'insulating';

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% time steps
Nt = 150;
dt = 60;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
temp_diff = max(kdiff.T(:)) - min(kdiff.T(:));

% check for test pass
if (temp_diff > insulating_heterog_diff)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end

% =========================================================================
% INSULATING - HETEROGENEOUS TEMP, PERFUSION
% =========================================================================

% set boundary condition
medium_perfused.boundary_condition = 'insulating';

% define arterial temperature
medium_perfused.blood_ambient_temperature = 0.5;  % [degC]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium_perfused, source, [], input_args{:});

% time steps
Nt = 100;
dt = 10;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
temp_diff = max(kdiff.T(:)) - min(kdiff.T(:));

% check for test pass
if (temp_diff > insulating_heterog_diff)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% create grid
Nx = 63;
Ny = 63;
Nz = 63;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define initial temperature distribution
source.T0 = 25 * makeBall(Nx, Ny, Nz, 0, 0, 0, 5);
source.T0 = smooth(source.T0, true);

% =========================================================================
% CONDUCTING - HETEROGENEOUS TEMP
% =========================================================================

% set boundary condition
medium.boundary_condition = 'conducting';

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% time steps
Nt = 150;
dt = 60;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
max_temp = abs(max(kdiff.T(:)));

% check for test pass
if (max_temp > conducting_max_temp)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end

% =========================================================================
% CONDUCTING - HETEROGENEOUS TEMP, PERFUSION
% =========================================================================

% set boundary condition
medium_perfused.boundary_condition = 'conducting';

% define arterial temperature
medium_perfused.blood_ambient_temperature = 0;  % [degC]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium_perfused, source, [], input_args{:});

% time steps
Nt = 100;
dt = 10;

% take time steps
kdiff.takeTimeStep(Nt, dt);

% check value has equalised
max_temp = abs(max(kdiff.T(:)));

% check for test pass
if (max_temp > conducting_max_temp)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    imagesc(squeeze(kdiff.T(:, :, round(Nz/2))));
    colorbar;
    axis image;
    subplot(2, 1, 2);
    plot(squeeze(kdiff.T(:, round(Ny/2), round(Nz/2))));
end