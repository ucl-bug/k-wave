function test_pass = kspaceFirstOrder3D_check_source_scaling_ux(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check the k-space source correction in
%     kspaceFirstOrder3D for a velocity source.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 13th December 2018
%     last update - 13th December 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set pass threshold
comparison_threshold = 5e-4;

% set pass variable
test_pass = true;

% =========================================================================
% CALCULATE GRID PARAMETERS
% =========================================================================

% create the computational grid
pml_size = 20;
Nx = 64;
Ny = 64;
Nz = 64;
dx = 50e-3/Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define the medium properties
medium.sound_speed = 1500;
medium.density = 1000;

% define a source mask
u_mask = zeros(Nx, Ny, Nz);
u_mask(1 + pml_size, 1 + pml_size:Ny - pml_size, 1 + pml_size:Nz - pml_size) = 1;

% define source frequency based on ppw
ppw = 3;
freq = medium.sound_speed / (ppw * dx);

% define ppp based on cfl, and force to be an integer
cfl = 0.5;
ppp = round(ppw / cfl);
dt = 1 / (ppp * freq);

% create the time array using an integer number of periods
t_end = 50e-6;
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% define time varying source signal
source_u = createCWSignals(kgrid.t_array, freq, 1, 0) / (medium.sound_speed * medium.density);

% define the different source conditions
source_k_correction.u_mask = u_mask;
source_k_correction.ux = source_u;

source_w_correction.u_mask = u_mask;
source_w_correction.ux = source_u;
source_w_correction.u_frequency_ref = freq;

source_no_correction.u_mask = u_mask;
source_no_correction.ux = source_u;
source_no_correction.u_mode = 'additive-no-correction';

% define sensor 
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(Nx - pml_size, Ny/2, Nz/2) = 1;

% only record the final few periods when the field is in steady state
record_periods = 3;
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

% =========================================================================
% RUN SIMULATIONS
% =========================================================================

% run the simulation
sensor_data_k_correction = kspaceFirstOrder3D(kgrid, medium, source_k_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size, ...
    'DataCast', 'single');

% get the amplitude
amp_k_correction = extractAmpPhase(sensor_data_k_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% run the simulation
sensor_data_w_correction = kspaceFirstOrder3D(kgrid, medium, source_w_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size, ...
    'DataCast', 'single');

% get the amplitude
amp_w_correction = extractAmpPhase(sensor_data_w_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% run the simulation
sensor_data_no_correction = kspaceFirstOrder3D(kgrid, medium, source_no_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size, ...
    'DataCast', 'single');

% get the amplitude
amp_no_correction = extractAmpPhase(sensor_data_no_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% reduce the cfl
fac = 20;
dt_fine = dt / fac;
Nt_fine = Nt * fac;
ppp_fine = ppp * fac;
kgrid.setTime(Nt_fine, dt_fine);

% recreate the source
source_no_correction.ux = createCWSignals(kgrid.t_array, freq, 1, 0) / (medium.sound_speed * medium.density);

% only record the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * ppp_fine + 1;

% run the simulation
sensor_data_fine = kspaceFirstOrder3D(kgrid, medium, source_no_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size, ...
    'DataCast', 'single');

% get the amplitude
amp_fine = extractAmpPhase(sensor_data_fine, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% =========================================================================
% COMPUTE ERROR
% =========================================================================

% compute errors
err_k_correction = abs(amp_k_correction - amp_fine) / amp_fine;
err_w_correction = abs(amp_w_correction - amp_fine) / amp_fine;
err_no_correction = abs(amp_no_correction - amp_fine) / amp_fine;

% check for test pass
if (err_k_correction > comparison_threshold) || (err_w_correction > comparison_threshold)
    test_pass = false;
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
   
    % create time axis
    t_array_fine = (0:(length(sensor_data_fine) - 1)) * dt_fine + dt_fine;
    t_array_coarse = (0:(length(sensor_data_k_correction) - 1)) * dt + dt;
    
    % plot
    figure;
    plot(t_array_fine, sensor_data_fine);
    hold on;    
    plot(t_array_coarse, sensor_data_k_correction, '.');
    plot(t_array_coarse, sensor_data_w_correction, 'o');
    plot(t_array_coarse, sensor_data_no_correction, 'x');
    legend('reference', 'k-space correction', 'w correction', 'no correction');
    title(['Errors = ' num2str(err_k_correction) ' (k), ' num2str(err_w_correction) ' (w), ' num2str(err_no_correction) ' (no)']);
    
end