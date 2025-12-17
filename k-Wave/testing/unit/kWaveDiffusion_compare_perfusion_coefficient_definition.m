function test_pass = kWaveDiffusion_compare_perfusion_coefficient_definition(~, plot_simulations) 
% DESCRIPTION:
%       Unit test to compare simulations defining the individual thermal
%       coefficients and the perfusion coefficient.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 18th August 2015
%       last update - 26th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015- Bradley Treeby

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
medium2 = medium1;

% define blood properties
medium1.blood_density        = 1060;	% [kg/m^3]
medium1.blood_specific_heat  = 3617;	% [J/(kg.K)]
medium1.blood_perfusion_rate = 0.01;  	% [1/s]
medium1.blood_ambient_temperature = 37;  % [degC]

% define equivalent perfusion coefficient
medium2.perfusion_coeff = medium1.blood_density .* medium1.blood_perfusion_rate .* medium1.blood_specific_heat ./ (medium1.density .* medium1.specific_heat);
medium2.blood_ambient_temperature = 37;  % [degC]

% set input args
input_args = {'PlotScale', [37, 38], 'PlotSim', plot_simulations};

% set initial temperature distribution
source.T0 = zeros(Nx, 1);
source.T0(Nx/2) = 1;
source.T0 = 37 + smooth(source.T0, true);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium1, source, [], input_args{:});

% take time steps
Nt = 50;
dt = 0.5;
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
Nt = 50;
dt = 0.5;
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

% =========================================================================
% 1D HETEROGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% clear old parameters
clear medium1 medium2

% define individual parameters
medium1.density              = 1079;     % [kg/m^3]
medium1.specific_heat        = 3540;     % [J/(kg.K)]
medium1.thermal_conductivity = 0.52;     % [W/(m.K)]
medium2 = medium1;

% define blood properties
blood_density_1        = 1060;  % [kg/m^3]
blood_density_2        = 1100;  % [kg/m^3]
blood_specific_heat_1  = 3617;  % [J/(kg.K)]
blood_specific_heat_2  = 3800;  % [J/(kg.K)]
blood_perfusion_rate_1 = 0.005; % [1/s]
blood_perfusion_rate_2 = 0.02;  % [1/s]
blood_temperature_1    = 30;
blood_temperature_2    = 37;

% define heterogeneous propeties
blood_density = blood_density_1 * ones(Nx, 1);
blood_density(kgrid.x > 0) = blood_density_2;

blood_specific_heat = blood_specific_heat_1 * ones(Nx, 1);
blood_specific_heat(kgrid.x > 0) = blood_specific_heat_2;

blood_perfusion_rate = blood_perfusion_rate_1 * ones(Nx, 1);
blood_perfusion_rate(kgrid.x > 0) = blood_perfusion_rate_2;

blood_temperature = blood_temperature_1 * ones(Nx, 1);
blood_temperature(kgrid.x > 0) = blood_temperature_2;

% define medium 1
medium1.blood_density        = blood_density;
medium1.blood_specific_heat  = blood_specific_heat;
medium1.blood_perfusion_rate = blood_perfusion_rate;
medium1.blood_ambient_temperature = blood_temperature;

% define equivalent perfusion coefficient
medium2.perfusion_coeff = blood_density .* blood_perfusion_rate .* blood_specific_heat ./ (medium1.density .* medium1.specific_heat);
medium2.blood_ambient_temperature = blood_temperature;  % [degC]

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
Nt = 50;
dt = 0.5;
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

% =========================================================================
% 2D HETEROGENEOUS SIMULATION
% =========================================================================

% create grid
Nx = 128;
Ny = 128;
dx = 1e-3;
dy = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define heterogeneous propeties
blood_density = blood_density_1 * ones(Nx, Ny);
blood_density(kgrid.y > 0) = blood_density_2;

blood_specific_heat = blood_specific_heat_1 * ones(Nx, Ny);
blood_specific_heat(kgrid.y > 0) = blood_specific_heat_2;

blood_perfusion_rate = blood_perfusion_rate_1 * ones(Nx, Ny);
blood_perfusion_rate(kgrid.y > 0) = blood_perfusion_rate_2;

blood_temperature = blood_temperature_1 * ones(Nx, Ny);
blood_temperature(kgrid.y > 0) = blood_temperature_2;

% define medium 1
medium1.blood_density        = blood_density;
medium1.blood_specific_heat  = blood_specific_heat;
medium1.blood_perfusion_rate = blood_perfusion_rate;
medium1.blood_ambient_temperature = blood_temperature;

% define equivalent perfusion coefficient
medium2.perfusion_coeff = blood_density .* blood_perfusion_rate .* blood_specific_heat ./ (medium1.density .* medium1.specific_heat);
medium2.blood_ambient_temperature = blood_temperature;

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

% define heterogeneous propeties
blood_density = blood_density_1 * ones(Nx, Ny, Nz);
blood_density(kgrid.y > 0) = blood_density_2;

blood_specific_heat = blood_specific_heat_1 * ones(Nx, Ny, Nz);
blood_specific_heat(kgrid.y > 0) = blood_specific_heat_2;

blood_perfusion_rate = blood_perfusion_rate_1 * ones(Nx, Ny, Nz);
blood_perfusion_rate(kgrid.y > 0) = blood_perfusion_rate_2;

blood_temperature = blood_temperature_1 * ones(Nx, Ny, Nz);
blood_temperature(kgrid.y > 0) = blood_temperature_2;

% define medium 1
medium1.blood_density        = blood_density;
medium1.blood_specific_heat  = blood_specific_heat;
medium1.blood_perfusion_rate = blood_perfusion_rate;
medium1.blood_ambient_temperature = blood_temperature;

% define equivalent perfusion coefficient
medium2.perfusion_coeff = blood_density .* blood_perfusion_rate .* blood_specific_heat ./ (medium1.density .* medium1.specific_heat);
medium2.blood_ambient_temperature = blood_temperature;

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