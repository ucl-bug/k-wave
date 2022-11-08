function test_pass = kWaveDiffusion_compare_diffusion_coefficient_definition(~, plot_simulations)
% DESCRIPTION:
%       Unit test to compare simulations defining the individual thermal
%       coefficients and the diffusion coefficient.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 18th August 2015
%       last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2019 Bradley Treeby

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
COMPARISON_THRESH   = 1e-13;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_simulations = true;
end

% =========================================================================
% 1D HOMOGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% define individual medium properties
medium1.density              = 1079;     % [kg/m^3]
medium1.thermal_conductivity = 0.52;     % [W/(m.K)]
medium1.specific_heat        = 3540;     % [J/(kg.K)]

% define equivalent diffusion coefficient
medium2.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);

% add density and specific heat for source scaling
medium3.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);
medium3.density              = medium1.density;
medium3.specific_heat        = medium1.specific_heat;

% set input args
input_args = {'PlotScale', [37, 38], 'PlotSim', plot_simulations};

% set initial temperature distribution
source.T0 = zeros(Nx, 1);
source.T0(Nx/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% number of time steps
Nt = 50;
dt = 0.5;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, 1);
source.Q(Nx/4) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% 2D HOMOGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% set initial temperature distribution
clear source;
source.T0 = zeros(Nx, Ny);
source.T0(Nx/2, Ny/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, Ny);
source.Q(Nx/4, Ny/2) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% 3D HOMOGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 64;
Ny = 64;
Nz = 64;
dx = 1e-3;
dy = 1e-3;
dz = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% set initial temperature distribution
clear source;
source.T0 = zeros(Nx, Ny, Nz);
source.T0(Nx/2, Ny/2, Nz/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, Ny, Nz);
source.Q(Nx/4, Ny/2, Nz/2) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% 1D HETEROGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% define heterogeneous medium properties
thermal_conducitivity_1 = 0.5;
thermal_conducitivity_2 = 0.2;

% clear old parameters
clear medium1 medium2

% define individual parameters
medium1.density              = 1079;     % [kg/m^3]
medium1.specific_heat        = 3540;     % [J/(kg.K)]
medium1.thermal_conductivity = thermal_conducitivity_1 .* ones(Nx, 1);
medium1.thermal_conductivity(kgrid.x > 0) = thermal_conducitivity_2;

% define equivalent diffusion coefficient
medium2.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);

% add density and specific heat for source scaling
medium3.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);
medium3.density              = medium1.density;
medium3.specific_heat        = medium1.specific_heat;

% set initial temperature distribution
clear source;
source.T0 = zeros(Nx, 1);
source.T0(Nx/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, 1);
source.Q(Nx/4) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% 2D HETEROGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define heterogeneous medium properties
medium1.thermal_conductivity = thermal_conducitivity_1 .* ones(Nx, Ny);
medium1.thermal_conductivity(kgrid.y > 0) = thermal_conducitivity_2;

% define equivalent diffusion coefficient
medium2.diffusion_coeff = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);

% add density and specific heat for source scaling
medium3.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);
medium3.density              = medium1.density;
medium3.specific_heat        = medium1.specific_heat;

% set initial temperature distribution
clear source;
source.T0 = zeros(Nx, Ny);
source.T0(Nx/2, Ny/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, Ny);
source.Q(Nx/4, Ny/2) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% 3D HETEROGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 64;
Ny = 64;
Nz = 64;
dx = 1e-3;
dy = 1e-3;
dz = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define heterogeneous medium properties
medium1.thermal_conductivity = thermal_conducitivity_1 .* ones(Nx, Ny, Nz);
medium1.thermal_conductivity(kgrid.y > 0) = thermal_conducitivity_2;

% define equivalent diffusion coefficient
medium2.diffusion_coeff = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);

% add density and specific heat for source scaling
medium3.diffusion_coeff      = medium1.thermal_conductivity ./ (medium1.density .* medium1.specific_heat);
medium3.density              = medium1.density;
medium3.specific_heat        = medium1.specific_heat;

% set initial temperature distribution
clear source;
source.T0 = zeros(Nx, Ny, Nz);
source.T0(Nx/2, Ny/2, Nz/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium2, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end

% -------------------------------------------------------------------------

% add a heat source
source.Q = zeros(Nx, Ny, Nz);
source.Q(Nx/4, Ny/2, Nz/2) = 1e5;
source.Q = smooth(source.Q, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T1 = kdiff.T;

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium3, source, [], input_args{:});

% take time steps
kdiff.takeTimeStep(Nt, dt);
T2 = kdiff.T;

% compute the maximum error
L_inf = max(abs(T2(:) - T1(:)));

% compute pass
if (L_inf > COMPARISON_THRESH)
    test_pass = false;
end