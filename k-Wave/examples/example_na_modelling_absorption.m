% Modelling Power Law Absorption Example
%
% This example describes the characteristics of the absorption and
% dispersion encapsulated by the k-Wave simulation functions.
%
% For a more detailed discussion of the absorption model used in k-Wave,
% see Treeby, B. E. and Cox, B. T., "Modeling power law absorption and
% dispersion for acoustic propagation using the fractional Laplacian," J.
% Acoust. Soc. Am., vol. 127, no. 5, pp. 2741-2748, 2010.
%
% author: Bradley Treeby
% date: 4th February 2011
% last update: 20th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% modify this parameter to run the different examples
example_number = 1;
% 1: simulation using kspaceSecondOrder 
% 2: simulation using kspaceFirstOrder1D with default CFL
% 3: simulation using kspaceFirstOrder1D with CFL = 0.05

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 1024;                      % [grid points]
x = 12.8e-3;                    % [m]
dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 4e-6;
if example_number < 3
    kgrid.makeTime(medium.sound_speed, [], t_end);
else
    kgrid.makeTime(medium.sound_speed, 0.05, t_end);
end

% create a spatial delta pulse
source_pos = Nx/4;              % [grid points]
source.p0 = zeros(Nx, 1);
source.p0(source_pos) = 1;

% define the sensor positions
source_sensor_dist = 0.5e-3;    % [m]
sensor_sensor_dist = 1e-3;      % [m]
sensor_pos_1 = source_pos + round(source_sensor_dist / dx);
sensor_pos_2 = source_pos + round((source_sensor_dist + sensor_sensor_dist) / dx);

% calculate discrete distance between the sensor positions
d = (sensor_pos_2 - sensor_pos_1) * dx;     % [m]
d_cm = d * 100;                             % [cm]

% index where the relative dispersion is defined
f_index = 30;

% create a Binary sensor mask
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos_1) = 1;
sensor.mask(sensor_pos_2) = 1;

% preallocate the storage variables
attenuation = zeros(3,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp_kk = zeros(3, floor(length(kgrid.t_array) / 2) + 1);

for loop = 1:3

    % define the absorption properties of the propagation medium
    switch loop
        case 1
            medium.alpha_coeff = 0.5;
            medium.alpha_power = 1.1;
        case 2         
            medium.alpha_coeff = 0.25;
            medium.alpha_power = 1.5;
        case 3           
            medium.alpha_coeff = 0.1;
            medium.alpha_power = 1.9;
    end

    % run the simulation without visualisation
    if example_number == 1
        sensor_data = kspaceSecondOrder(kgrid, medium, source, sensor, 'PlotSim', false);   
    else
        sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, 'PlotSim', false);
    end
    
    % calculate the amplitude and phase spectrum at the two sensor
    % positions
    [~, as1, ps1] = spect(sensor_data(1, :), 1/kgrid.dt);
    [f, as2, ps2] = spect(sensor_data(2, :), 1/kgrid.dt);
    
    % calculate the attenuation from the amplitude spectrums
    attenuation(loop, :) = -20 * log10(as2 ./ as1) ./ d_cm;
    
    % calculate the corresponding theoretical attenuation in dB/cm
    attenuation_th(loop, :) = medium.alpha_coeff .* (f .* 1e-6).^medium.alpha_power;
    
    % calculate the dispersion (dependence of the sound speed on frequency)
    % from the phase spectrums
    cp(loop, :) = 2 .* pi .* f .* d ./ (unwrap(ps1) - unwrap(ps2));
    
    % calculate the corresponding theoretical dispersion using the
    % Kramers-Kronig relation for power law absorption
    cp_kk(loop, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(loop, f_index), db2neper(medium.alpha_coeff, medium.alpha_power), medium.alpha_power);
    
end

% =========================================================================
% VISUALISATION
% =========================================================================

% set downsampling factor so there is sufficient space between plot markers
ds = 8;

% plot the attenuation
figure;
f_max = 50;     % [MHz]
plot(f(1:ds:end) .* 1e-6, attenuation(:, 1:ds:end), 'ko', f .* 1e-6, attenuation_th, 'k-');
set(gca, 'XLim', [0, f_max]);
box on;
xlabel('Frequency [MHz]');
ylabel('\alpha [dB/cm]');

% label the plots
text(40, 160, 'y = 1.9');
text(40, 90, 'y = 1.5');
text(40, 40, 'y = 1.1');

% plot the dispersion
figure
plot(f(1:ds:end) .* 1e-6, cp(:, 1:ds:end), 'ko', f .* 1e-6, cp_kk, 'k-');
set(gca, 'XLim', [0, f_max]);
box on;
xlabel('Frequency [MHz]');
ylabel('C_p [m/s]');

% label the plots
text(40, 1517, 'y = 1.1');
text(40, 1509, 'y = 1.5');
text(40, 1500, 'y = 1.9');