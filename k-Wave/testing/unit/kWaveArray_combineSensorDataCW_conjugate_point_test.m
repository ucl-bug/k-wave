function test_pass = kWaveArray_combineSensorDataCW_conjugate_point_test(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to verify the combineSensorDataCW function for hologram 
%       elements by testing the acoustic reciprocity principle for a point
%       source and hologram with a single point or multiple points.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 25th March 2025
%       last update - 25th March 2025
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2025- Bradley Treeby

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

% set pass variable
test_pass = true;

% set comparison threshold (percentage error)
comparison_thresh = 1e-3;

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% compute points per temporal period
source_freq = 250e3;  % [Hz]
ppw = 4;
c0 = 1500;  % [m/s]
cfl = 0.3;
dx = c0 / (ppw * source_freq);   % [m]
PPP = round(ppw / cfl);

% grid properties
Nx = 64;
Ny = 32;
Nz = 32;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% medium properties
medium.sound_speed = c0;

% create the time array using an integer number of points per period
dt = 1 / (PPP * source_freq);
t_end = Nx * dx / c0;
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

for test_ind = 1:2

    % =====================================================================
    % FORWARD SIMULATION: POINT SOURCE AND HOLOGRAM SENSOR
    % =====================================================================
    
    % create point source 
    source_pt.p_mask = zeros(Nx, Ny, Nz);
    source_pt.p_mask(Nx/2 + 1, Ny/2 + 1, Nz/2 + 1) = 1;
    
    % create CW source
    source_phase = rand * pi;
    source_pt.p = createCWSignals(kgrid.t_array, source_freq, 1, source_phase);
    
    % create hologram array for sensing
    karray = kWaveArray();
    
    % add hologram element to the array
    hologram_position = [-kgrid.x_size/4, 0, 0].';
    hologram_area = (kgrid.dx)^2;
    switch test_ind
        case 1
            hologram_phase = rand * pi;
            karray.addHologramElement(hologram_position, hologram_position, ...
                source_freq, 1, hologram_phase, hologram_area);
        case 2
            hologram_integration_points = [-kgrid.x_size/4, -dx, 0; -kgrid.x_size/4, 0, 0; -kgrid.x_size/4, dx, 0].';
            hologram_amp = rand(1, 3);
            hologram_phase = rand(1, 3) * pi;
            karray.addHologramElement(hologram_position, hologram_integration_points, ...
                source_freq, hologram_amp, hologram_phase, hologram_area);
    end
    
    % set sensor mask using hologram array
    sensor_hologram.mask = karray.getArrayBinaryMask(kgrid);
    
    % configure sensor to record time-varying pressure
    sensor_hologram.record = {'p'};
    
    % run first simulation with point source and hologram sensor
    sensor_data_pt_src = kspaceFirstOrder3D(kgrid, medium, source_pt, sensor_hologram, ...
        'PlotSim', plot_simulations, ...
        'PlotScale', 'auto');
    
    % extract amplitude from the sensor data
    [amp, phase] = extractAmpPhase(sensor_data_pt_src.p(:, end - PPP + 1:end), 1/kgrid.dt, source_freq, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    
    % as we want to reference against the source phase, account for phase
    % asymmetry between createCWSignals (uses sin) and extractAmpPhase
    % (implicitly assumes cos)
    phase = phase - pi/2;
    
    % use combineSensorDataCW to get complex amplitude at the array element
    sensor_data_complex = karray.combineSensorDataCW(kgrid, amp .* exp(1i * phase));
    
    % =====================================================================
    % REVERSE SIMULATION: HOLOGRAM SOURCE AND POINT SENSOR 
    % =====================================================================
    
    % set source mask using hologram array
    source_hologram.p_mask = karray.getArrayBinaryMask(kgrid);
    
    % get distributed source signal for CW input
    source_hologram.p = karray.getDistributedSourceSignalCW(kgrid, abs(sensor_data_complex), -angle(sensor_data_complex));
    
    % create point sensor at the original source location
    sensor_pt.mask = source_pt.p_mask;
    sensor_pt.record = {'p'};
    
    % run second simulation with hologram source and point sensor
    sensor_data_holo_src = kspaceFirstOrder3D(kgrid, medium, source_hologram, sensor_pt, ...
        'PlotSim', plot_simulations, ...
        'PlotScale', 'auto');
    
    % extract amplitude from the sensor data, accounting for sin/cos
    [~, phase_pt] = extractAmpPhase(sensor_data_holo_src.p(:, end - PPP + 1:end), 1/kgrid.dt, source_freq, ...
        'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    phase_pt = phase_pt - pi/2;
    
    % conjugate phase a final time to match the input phase
    phase_pt = -phase_pt;
    
    % =========================================================================
    % COMPARE RESULTS AND CALCULATE ERROR
    % =========================================================================
    
    % due to reciprocity, the phase should match the original source phase
    reciprocity_error_phase = abs(mod(phase_pt - source_phase, 2*pi)) * 180/pi;
    reciprocity_error_phase = min(reciprocity_error_phase, 360 - reciprocity_error_phase);
    fprintf('Phase error: %.2f degrees\n', reciprocity_error_phase);
    
    % check if error is within threshold
    if (reciprocity_error_phase > comparison_thresh)
        test_pass = false;
        fprintf('TEST FAILED: Error exceeds threshold\n');
    else
        fprintf('TEST PASSED: Error within threshold\n');
    end
    
    if plot_comparisons
        figure;
        subplot(1, 2, 1);
        plot(kgrid.t_array, source_pt.p ./ max(source_pt.p));
        hold on;
        plot(kgrid.t_array, sensor_data_holo_src.p ./ max(sensor_data_holo_src.p), '--');
        legend('Input Forward', '- Output Back');
    
        subplot(1, 2, 2);
        plot(source_pt.p(:, end - PPP + 1:end) ./ max(source_pt.p));
        hold on;
        plot(sensor_data_holo_src.p(:, end - PPP + 1:end) ./ max(sensor_data_holo_src.p), '--');
        legend('Input Forward', 'Output Back');
    
        scaleFig(2, 1);
    end

end