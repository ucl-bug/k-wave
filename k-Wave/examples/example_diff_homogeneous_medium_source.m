% Constant Rate Of Heat Deposition Example
%
% This example provides a simple demonstration of using k-Wave to model the
% diffusion of heat within a two-dimensional homogeneous medium with a
% constant source term (volume rate of heat deposition). It builds on the
% Heat Diffusion In A Homogeneous Medium Example.
%
% author: Bradley Treeby
% date: 15th January 2017
% last update: 28th April 2017
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

clearvars;

% =========================================================================
% SIMULATION USING KWAVEDIFFUSION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 1e-3;          % grid point spacing in the x direction [m]
dy = 1e-3;          % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define medium properties related to diffusion
medium.density                      = 1079;     % [kg/m^3]
medium.thermal_conductivity         = 0.52;     % [W/(m.K)]
medium.specific_heat                = 3540;     % [J/(kg.K)]

% define medium properties related to perfusion
medium.blood_density                = 1060;     % [kg/m^3]
medium.blood_specific_heat          = 3617;     % [J/(kg.K)]
medium.blood_perfusion_rate         = 0.01;     % [1/s]
medium.blood_ambient_temperature    = 37;       % [degC]

% set initial temperature distribution to be constant [degC]
source.T0 = 37;

% set Gaussian volume rate of heat deposition [W/m^3]
width = 4*dx;
source.Q = 2e6 .* exp( -(kgrid.x ./ width).^2 - (kgrid.y ./ width).^2 );

% set input args
input_args = {'PlotScale', [37, 50]};

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, [], input_args{:});

% take time steps (temperature can be accessed as kdiff.T)
Nt = 300;
dt = 0.5;
kdiff.takeTimeStep(Nt, dt);

% plot the current temperature field
figure;
kdiff.plotTemp;

% =========================================================================
% SIMULATION USING BIOHEATEXACT
% =========================================================================

% calculate perfusion coefficient from the medium parameters
P = medium.blood_density .* medium.blood_perfusion_rate .* ...
    medium.blood_specific_heat ./ (medium.density .* medium.specific_heat);

% calculate diffusivity from the medium parameters
D = medium.thermal_conductivity / (medium.density * medium.specific_heat);

% calculate normalised heat source
S = source.Q ./ (medium.density .* medium.specific_heat);

% compute Green's function solution using bioheatExact
T_exact = bioheatExact(source.T0 .* ones(Nx, Ny), S, ...
    [D, P, medium.blood_ambient_temperature], kgrid.dx, Nt * dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the thermal dose and lesion map
figure;

subplot(1, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.cem43, [0, 240]);
h = colorbar;
xlabel(h, '[CEM43]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Thermal Dose');

subplot(1, 2, 2);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.lesion_map, [0, 1]);
colorbar;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Ablated Tissue');

colormap(flipud(hot(256)));
scaleFig(1.5, 1);

% plot final temperature and thermal dose
figure;

subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.T);
h = colorbar;
xlabel(h, '[^\circC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Final Temperatue (kWaveDiffusion)');

subplot(2, 2, 2);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, T_exact);
h = colorbar;
xlabel(h, '[^\circC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Final Temperature (bioheatExact)');

subplot(2, 2, 3);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, abs(T_exact - kdiff.T));
colorbar;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Difference');

colormap(jet(256));