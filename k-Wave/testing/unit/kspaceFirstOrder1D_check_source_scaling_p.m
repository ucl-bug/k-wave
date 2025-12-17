function test_pass = kspaceFirstOrder1D_check_source_scaling_p(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check the k-space source correction in
%     kspaceFirstOrder1D for a pressure source.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th May 2018
%     last update - 13th December 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018- Bradley Treeby

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
comparison_threshold = 1e-3;

% set pass variable
test_pass = true;

% =========================================================================
% CALCULATE GRID PARAMETERS
% =========================================================================

% create the computational grid
pml_size = 20;
Nx = 64;
dx = 50e-3/Nx;
kgrid = kWaveGrid(Nx, dx);

% define the medium properties
medium.sound_speed = 1500;

% define a source mask
p_mask = zeros(Nx, 1);
p_mask(pml_size + 1) = 1;

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
source_p = createCWSignals(kgrid.t_array, freq, 1, 0);

% define the different source conditions
source_k_correction.p_mask = p_mask;
source_k_correction.p = source_p;

source_w_correction.p_mask = p_mask;
source_w_correction.p = source_p;
source_w_correction.p_frequency_ref = freq;

source_no_correction.p_mask = p_mask;
source_no_correction.p = source_p;
source_no_correction.p_mode = 'additive-no-correction';

% define sensor 
sensor.mask = zeros(Nx, 1);
sensor.mask(Nx - pml_size) = 1;

% only record the final few periods when the field is in steady state
record_periods = 3;
sensor.record_start_index = kgrid.Nt - record_periods * ppp + 1;

% =========================================================================
% RUN SIMULATIONS
% =========================================================================

% run the simulation
sensor_data_k_correction = kspaceFirstOrder1D(kgrid, medium, source_k_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size);

% get the amplitude
amp_k_correction = extractAmpPhase(sensor_data_k_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% run the simulation
sensor_data_w_correction = kspaceFirstOrder1D(kgrid, medium, source_w_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size);

% get the amplitude
amp_w_correction = extractAmpPhase(sensor_data_w_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% run the simulation
sensor_data_no_correction = kspaceFirstOrder1D(kgrid, medium, source_no_correction, sensor, ...
    'PlotSim', plot_simulations, ...
    'PMLSize', pml_size);

% get the amplitude
amp_no_correction = extractAmpPhase(sensor_data_no_correction, 1/kgrid.dt, freq, 'Window', 'Rectangular', 'FFTPadding', 1);

% =========================================================================
% COMPUTE ERROR
% =========================================================================

% compute errors (amplitude should be 1)
err_k_correction = amp_k_correction - 1;
err_w_correction = amp_w_correction - 1;
err_no_correction = amp_no_correction - 1;

% check for test pass
if (err_k_correction > comparison_threshold) || (err_w_correction > comparison_threshold)
    test_pass = false;
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
   
    figure;
    plot(sensor_data_k_correction);
    hold on;
    plot(sensor_data_w_correction, '--');
    plot(sensor_data_no_correction, ':');
    legend('k-space correction', 'w correction', 'no correction');
    title(['Errors = ' num2str(err_k_correction) ' (k), ' num2str(err_w_correction) ' (w), ' num2str(err_no_correction) ' (no)']);
    
end