function test_pass = kWaveDiffusion_compare_cem_with_exact(~, plot_simulations) 
% DESCRIPTION:
%       Unit test to check the cumulative equivalent minutes calculated by
%       kWaveDiffusion is correct.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 18th August 2015
%       last update - 11th December 2017
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

% set pass variable
test_pass = true;

% set comparison threshold
COMPARISON_THRESH   = 1e-14;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_simulations = true;
end

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

% set constant initial temperature distribution
set_temp = 50;
source.T0 = set_temp;

% set input args
input_args = {'PlotScale', [37, 38], 'PlotSim', plot_simulations};

% Thermal dose formula for constant temperature
%   cem43 = t * R^(43 - T)
%   R = 0.5  for  T >= 43
%   R = 0.25 for  43 > T >= 39

% calculate time to give a thermal dose of 200 minutes at 50 degrees
t200 = (200*60) / 0.5^(43 - set_temp); % [s]

% calculate time to give a thermal dose of 300 minutes at 50 degrees
t300 = (300*60) / 0.5^(43 - set_temp); % [s]

% compute equivalent number of time steps
dt = 0.5;
Nt200 = round(t200/dt);
Nt300 = round(t300/dt);

% calculate exact thermal dose
ref_cem_1 = Nt200 * dt / 60 * 0.5^(43 - set_temp); % [mins]
ref_cem_2 = Nt300 * dt / 60 * 0.5^(43 - set_temp); % [mins]

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt200, dt);

% extract cem and lesion volume
cem_1 = kdiff.cem43(1);
vol_1 = kdiff.lesion_size;

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt300, dt);

% extract cem and lesion volume
cem_2 = kdiff.cem43(1);
vol_2 = kdiff.lesion_size;

% compute pass
ref_vol_1 = 0;
ref_vol_2 = Nx*dx;
if (cem_1 - ref_cem_1) > COMPARISON_THRESH
    test_pass = false;
end
if (cem_2 - ref_cem_2) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_1 - ref_vol_1) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_2 - ref_vol_2) > COMPARISON_THRESH
    test_pass = false;
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

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt200, dt);

% extract cem and lesion volume
cem_1 = kdiff.cem43(1);
vol_1 = kdiff.lesion_size;

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt300, dt);

% extract cem and lesion volume
cem_2 = kdiff.cem43(1);
vol_2 = kdiff.lesion_size;

% compute pass
ref_vol_1 = 0;
ref_vol_2 = Nx*dx*Ny*dy;
if (cem_1 - ref_cem_1) > COMPARISON_THRESH
    test_pass = false;
end
if (cem_2 - ref_cem_2) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_1 - ref_vol_1) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_2 - ref_vol_2) > COMPARISON_THRESH
    test_pass = false;
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

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt200, dt);

% extract cem and lesion volume
cem_1 = kdiff.cem43(1);
vol_1 = kdiff.lesion_size;

% create kWaveDiffusion object:
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt300, dt);

% extract cem and lesion volume
cem_2 = kdiff.cem43(1);
vol_2 = kdiff.lesion_size;

% compute pass
ref_vol_1 = 0;
ref_vol_2 = Nx*dx*Ny*dy*Nz*dz;
if (cem_1 - ref_cem_1) > COMPARISON_THRESH
    test_pass = false;
end
if (cem_2 - ref_cem_2) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_1 - ref_vol_1) > COMPARISON_THRESH
    test_pass = false;
end
if (vol_2 - ref_vol_2) > COMPARISON_THRESH
    test_pass = false;
end