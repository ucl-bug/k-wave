function test_pass = kspaceFirstOrderAS_check_source_scaling_ux(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check the k-space source correction in
%     kspaceFirstOrderAS for a velocity source.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 18th December 2018
%     last update - 18th December 2018
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
comparison_threshold_k = 2e-4;
comparison_threshold_w = 5e-2;

% set pass variable
test_pass = true;

% =========================================================================
% CALCULATE GRID PARAMETERS
% =========================================================================

% create grid
Nx = 54;
Ny = 134;
dx = 3.75e-4;
dy = dx;
kgrid_coarse = kWaveGrid(Nx, dx, Ny, dy);
kgrid_fine   = kWaveGrid(Nx, dx, Ny, dy);

% assign medium properties
medium.sound_speed = 1500;
medium.density = 1000;

% create time array
t_end  = 40e-6;
cfl = 0.5;
kgrid_coarse.makeTime(medium.sound_speed, cfl, t_end);

% create second time array
sc_fac = 10;
cfl = cfl / sc_fac;
kgrid_fine.makeTime(medium.sound_speed, cfl, t_end);

% tone burst
freq = 1e6;
source_sig_coarse = toneBurst(1/kgrid_coarse.dt, freq, 3) / (medium.sound_speed * medium.density);
source_sig_fine   = toneBurst(1/kgrid_fine.dt, freq, 3) / (medium.sound_speed * medium.density);

% time index where sensor data is split into direct and edge waves
split_index = 3200;

% create piston source mask
u_mask = zeros(Nx, Ny);
u_mask(1, 1:80) = 1;

% define the different source conditions
source_fine.u_mask = u_mask;
source_fine.ux = source_sig_fine;

source_k_correction.u_mask = u_mask;
source_k_correction.ux = source_sig_coarse;

source_w_correction.u_mask = u_mask;
source_w_correction.ux = source_sig_coarse;
source_w_correction.u_frequency_ref = freq;

source_no_correction.u_mask = u_mask;
source_no_correction.ux = source_sig_coarse;
source_no_correction.u_mode = 'additive-no-correction';

% set sensor mask
sensor.mask = zeros(Nx, Ny);
sensor.mask(end, 1) = 1;

% set input options
input_args = {...
    'PlotSim', plot_simulations, ...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    };

% =========================================================================
% RUN SIMULATIONS
% =========================================================================

% set symmetries to test
if exist('dtt1D', 'file') 
    symmetry_array = {'WSWA-FFT', 'WSWS-FFT', 'WSWA', 'WSWS'};
else
    symmetry_array = {'WSWA-FFT', 'WSWS-FFT'};
end

% loop through symmetries
for ind = 1:length(symmetry_array)

    % set the symmetry
    symmetry = symmetry_array{ind};
    
    % run the simulations
    sensor_data_k = kspaceFirstOrderAS(kgrid_coarse, medium, source_k_correction, sensor, ...
        input_args{:}, ...
        'RadialSymmetry', symmetry);

    sensor_data_w = kspaceFirstOrderAS(kgrid_coarse, medium, source_w_correction, sensor, ...
        input_args{:}, ...
        'RadialSymmetry', symmetry);
    
    sensor_data_no_corr = kspaceFirstOrderAS(kgrid_coarse, medium, source_no_correction, sensor, ...
        input_args{:}, ...
        'RadialSymmetry', symmetry);
    
    sensor_data_ref = kspaceFirstOrderAS(kgrid_fine, medium, source_fine, sensor, ...
        input_args{:}, ...
        'RadialSymmetry', symmetry);
    
    % upsample simulations
    sensor_data_k_us        = interpft(sensor_data_k,       length(sensor_data_k)       * sc_fac * 2);
    sensor_data_w_us        = interpft(sensor_data_w,       length(sensor_data_w)       * sc_fac * 2);
    sensor_data_no_corr_us  = interpft(sensor_data_no_corr, length(sensor_data_no_corr) * sc_fac * 2);
    sensor_data_ref_us      = interpft(sensor_data_ref,     length(sensor_data_ref)     * 2);

    % account for time staggering
    sensor_data_k_us        = circshift(sensor_data_k_us, 2 * sc_fac);
    sensor_data_w_us        = circshift(sensor_data_w_us, 2 * sc_fac);
    sensor_data_no_corr_us  = circshift(sensor_data_no_corr_us, 2 * sc_fac);
    sensor_data_ref_us      = circshift(sensor_data_ref_us, 2); 

    % trim to the same length
    sensor_data_k_us        = sensor_data_k_us(1:length(sensor_data_ref_us));
    sensor_data_w_us        = sensor_data_w_us(1:length(sensor_data_ref_us));
    sensor_data_no_corr_us  = sensor_data_no_corr_us(1:length(sensor_data_ref_us));

    % normalise
    norm_fac                = max(sensor_data_ref_us(:));
    sensor_data_ref_us      = sensor_data_ref_us ./ norm_fac;
    sensor_data_k_us        = sensor_data_k_us ./ norm_fac;
    sensor_data_w_us        = sensor_data_w_us ./ norm_fac;
    sensor_data_no_corr_us  = sensor_data_no_corr_us ./ norm_fac;

    % get error, separating into direct wave and edge wave
    err_k_direct            = max(abs(sensor_data_k_us(1:split_index) - sensor_data_ref_us(1:split_index)));
    err_w_direct            = max(abs(sensor_data_w_us(1:split_index) - sensor_data_ref_us(1:split_index)));
    err_no_corr_direct      = max(abs(sensor_data_no_corr_us(1:split_index) - sensor_data_ref_us(1:split_index)));
    
    err_k_edge              = max(abs(sensor_data_k_us(split_index + 1:end) - sensor_data_ref_us(split_index + 1:end)));
    err_w_edge              = max(abs(sensor_data_w_us(split_index + 1:end) - sensor_data_ref_us(split_index + 1:end)));
    err_no_corr_edge        = max(abs(sensor_data_no_corr_us(split_index + 1:end) - sensor_data_ref_us(split_index + 1:end)));
    
    % check for test pass
    if (err_k_direct > comparison_threshold_k) || ...
            (err_w_direct > comparison_threshold_w) || ...
            (err_k_direct > err_w_direct) || ...
            (err_w_direct > err_no_corr_direct) || ...
            (err_k_edge > err_w_edge) || ...
            (err_w_edge > err_no_corr_edge)
        test_pass = false;
    end
    
    if plot_comparisons
        
        % plot time series
        figure;
        plot(sensor_data_ref_us);
        hold on;
        plot(sensor_data_k_us, '.');
        plot(sensor_data_w_us, '.');
        plot(sensor_data_no_corr_us, '.');
        legend('reference', 'k-space correction', 'w correction', 'no correction');
        ylabel('sensor data');
        xlabel('time index');
        title({symmetry, ...
               ['Direct Wave Errors = ' num2str(err_k_direct) ' (k), ' num2str(err_w_direct) ' (w), ' num2str(err_no_corr_direct) ' (no)'], ...
               ['Edge Wave Errors = ' num2str(err_k_edge) ' (k), ' num2str(err_w_edge) ' (w), ' num2str(err_no_corr_edge) ' (no)']});
        
    end
    
end