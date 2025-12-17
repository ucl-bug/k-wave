function test_pass = kspaceFirstOrder_compare_2D_AS_cylindrical_waves(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare a cylindrical wave in 2D and using the
%     axisymmetric code to catch any coding bugs between the two codes. A
%     homogeneous medium is used.
%
%     8 tests are performed:
%
%         1.  linear    + lossless + source.p (additive)
%         2.  linear    + lossless + source.p (dirichlet)
%         3.  linear    + stokes   + source.p (additive)
%         4.  linear    + stokes   + source.p (dirichlet)
%         5.  nonlinear + lossless + source.p (additive)
%         6.  nonlinear + lossless + source.p (dirichlet)
%         7.  nonlinear + stokes   + source.p (additive)
%         8.  nonlinear + stokes   + source.p (dirichlet)
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 2nd November 2017
%     last update - 29th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017- Bradley Treeby

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

%#ok<*NOPRT>
%#ok<*UNRCH>

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set additional literals to give further permutations of the test
STAGGERED_GRID      = true;
USE_KSPACE          = true;
COMPARISON_THRESH   = 5e-3;
USE_PML             = true;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 32;
Ny = 128;
dx = 25e-3/Nx;    
dy = dx; 

% define PML properties
PML_size    = 20;
if USE_PML
    PML_alpha   = 2;
else
    PML_alpha   = 0; 
end
    
% define medium properties
c0      = 1500;
rho0    = 1000;
BonA0   = 10;
alpha0  = 5;

% define source signal
CFL     = 0.2;
PPW     = 6;
amp     = 1;
phase   = 0;

% other parameters
num_record_periods = 3;

% define optional inputs
input_args_2D = {'PMLSize', PML_size, 'PMLAlpha', PML_alpha, 'PMLInside', false, ...
    'UsekSpace', USE_KSPACE, 'UseSG', STAGGERED_GRID, ...
    'PlotSim', plot_simulations, 'PlotScale', 'auto', 'PlotPML', false};

input_args_as = {'PMLSize', [5, PML_size], 'PMLAlpha', [0, PML_alpha], ...
    'UsekSpace', USE_KSPACE, 'UseSG', STAGGERED_GRID, ...
    'PlotSim', plot_simulations, 'PlotScale', 'auto'};

% set pass variable
test_pass = true;

% test names
test_names = {...
    'linear + lossless + source.p (additive)', ...
    'linear + lossless + source.p (dirichlet)', ...
    'linear + stokes + source.p (additive)', ...
    'linear + stokes + source.p (dirichlet)', ...
    'nonlinear + lossless + source.p (additive)', ...
    'nonlinear + lossless + source.p (dirichlet)', ...
    'nonlinear + stokes + source.p (additive)', ...
    'nonlinear + stokes + source.p (dirichlet)', ...
    };

% lists used to set properties
dirichlet_tests = 2:2:8;

% =========================================================================
% SIMULATIONS
% =========================================================================

% create the computational grid    
kgrid   = kWaveGrid(Nx, dx, Ny, dy);

% calculate source frequency
freq    = c0 / (PPW * dx);

% define time array to give an integer number of points per period
dt      = CFL * dx / c0;
PPP     = round(1 / (dt * freq));
dt      = 1 / (freq * PPP);

% set time to propagate across grid
t_end   = 1.5 * sqrt( (Nx * dx).^2 + (Ny * dy).^2 ) / c0;
Nt      = round(t_end / dt);

% assign to kgrid
kgrid.setTime(Nt, dt);

% create source signal
source_sig = createCWSignals(kgrid.t_array, freq, amp, phase);

% define line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2, 1:end - PML_size) = 1;
sensor.record = {'p'};

% only record the last few periods in steady state
sensor.record_start_index = Nt - num_record_periods * PPP + 1;

% loop through tests
for test_num = 1:8

    % clear structures
    clear source medium
    
    % update command line
    disp(['Running Test: ' test_names{test_num}]);    
    
    % assign medium properties
    medium.sound_speed          = c0;
    medium.density              = rho0;
    switch test_num
        case {1, 2}
            
            % linear + lossless
            
        case {3, 4}
            
            % linear + stokes
            medium.alpha_coeff  = alpha0;
            medium.alpha_power	= 2;    
            
        case {5, 6}
            
            % nonlinear + lossless
            medium.BonA         = BonA0;
            
        case {7, 8}
            
            % nonlinear + stokes
            medium.alpha_coeff  = alpha0;
            medium.BonA         = BonA0;
            medium.alpha_power  = 2;
            
    end
    
    % ------------------------------------
    % AXISYMMETRIC SOURCE
    % ------------------------------------
    
    % define source
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(:, 1) = 1;
    source.p = source_sig;
    if any(dirichlet_tests == test_num)
        source.p_mode = 'dirichlet';
    end        
    
    % define source mask
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(:, 1) = 1;
   
    % ------------------------------------
    % AXISYMMETRIC SIMULATION: WSWS-FFT
    % ------------------------------------

    % run the simulation
    sensor_data_as_wsws_fft = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
        input_args_as{:}, 'RadialSymmetry', 'WSWS-FFT');
        
    % extract amplitude and phase from recorded data
    sensor_data_as_wsws_fft = extractAmpPhase(sensor_data_as_wsws_fft.p, 1/kgrid.dt, freq, 'Dim', 2);    
    
    % ------------------------------------
    % AXISYMMETRIC SIMULATION: WSWA-FFT
    % ------------------------------------

    % run the simulation
    sensor_data_as_wswa_fft = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
        input_args_as{:}, 'RadialSymmetry', 'WSWA-FFT');
        
    % extract amplitude and phase from recorded data
    sensor_data_as_wswa_fft = extractAmpPhase(sensor_data_as_wswa_fft.p, 1/kgrid.dt, freq, 'Dim', 2);
    
    % ------------------------------------
    % AXISYMMETRIC SIMULATION: WSWS
    % ------------------------------------

    if exist('dtt1D', 'file')
    
        % run the simulation
        sensor_data_as_wsws = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
            input_args_as{:}, 'RadialSymmetry', 'WSWS');

        % extract amplitude and phase from recorded data
        sensor_data_as_wsws = extractAmpPhase(sensor_data_as_wsws.p, 1/kgrid.dt, freq, 'Dim', 2);
        
    end
    
    % ------------------------------------
    % AXISYMMETRIC SIMULATION: WSWA
    % ------------------------------------

    if exist('dtt1D', 'file')
    
        % run the simulation
        sensor_data_as_wswa = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
            input_args_as{:}, 'RadialSymmetry', 'WSWA');

        % extract amplitude and phase from recorded data
        sensor_data_as_wswa = extractAmpPhase(sensor_data_as_wswa.p, 1/kgrid.dt, freq, 'Dim', 2);
        
    end    
    
    % ------------------------------------
    % 2D SIMULATION: X
    % ------------------------------------
    
    % define source
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(Nx/2, 1) = 1;
    source.p = source_sig;
    if any(dirichlet_tests == test_num)
        source.p_mode = 'dirichlet';
    end        
    
    % run the simulation
    sensor_data_2D_x = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        input_args_2D{:});

    % extract amplitude and phase from recorded data
    sensor_data_2D_x = extractAmpPhase(sensor_data_2D_x.p, 1/kgrid.dt, freq, 'Dim', 2);
    
    % ------------------------------------
    % COMPARISON
    % ------------------------------------

    % radial axis and normalisation point
    x = (0:length(sensor_data_as_wsws_fft) - 1) * dx;
    [~, ind] = findClosest(x, 10e-3);
    
    % normalise
    sensor_data_2D_x        = sensor_data_2D_x          ./ sensor_data_2D_x(ind);
    sensor_data_as_wsws_fft = sensor_data_as_wsws_fft   ./ sensor_data_as_wsws_fft(ind);
    sensor_data_as_wswa_fft = sensor_data_as_wswa_fft   ./ sensor_data_as_wswa_fft(ind);
    if exist('dtt1D', 'file')
        sensor_data_as_wsws = sensor_data_as_wsws       ./ sensor_data_as_wsws(ind);
        sensor_data_as_wswa = sensor_data_as_wswa       ./ sensor_data_as_wswa(ind);    
    end
    
    % compute errors
    ref_max = max(abs(sensor_data_2D_x(ind:end)));
    diff_2D_as_wsws_fft = max(abs(sensor_data_2D_x(ind:end) - sensor_data_as_wsws_fft(ind:end))) / ref_max 
    diff_2D_as_wswa_fft = max(abs(sensor_data_2D_x(ind:end) - sensor_data_as_wswa_fft(ind:end))) / ref_max 
    if exist('dtt1D', 'file')
        diff_2D_as_wsws = max(abs(sensor_data_2D_x(ind:end) - sensor_data_as_wsws(ind:end))) / ref_max 
        diff_2D_as_wswa = max(abs(sensor_data_2D_x(ind:end) - sensor_data_as_wswa(ind:end))) / ref_max
    else
        diff_2D_as_wsws = 0;
        diff_2D_as_wswa = 0;
        disp('WARNING: ''RadialSymmetry'' options ''WSWS'' and ''WSWA'' not tested, as dtt library is not present on path');
    end
    
    % check for test pass
    if (diff_2D_as_wsws_fft > COMPARISON_THRESH) || ...
       (diff_2D_as_wswa_fft > COMPARISON_THRESH) || ...
       (diff_2D_as_wsws > COMPARISON_THRESH) || ...
       (diff_2D_as_wswa > COMPARISON_THRESH)      
        test_pass = false;
    end    
    
    % ------------------------------------
    % PLOTTING
    % ------------------------------------    
    
    if plot_comparisons
        
        % fit 1/sqrt(r)
        decay = 1 ./ sqrt(x);
        decay = decay ./ decay(ind);
        
        figure;
        
        subplot(4, 1, 1)
        plot(x, sensor_data_2D_x);
        hold on;
        plot(x, sensor_data_as_wswa_fft, '--');
        plot(x, decay, ':')
        title(['AS WSWA-FFT, L_{inf} = ' num2str(diff_2D_as_wswa_fft)]);
        
        if exist('dtt1D', 'file')
            subplot(4, 1, 2)
            plot(x, sensor_data_2D_x);
            hold on;
            plot(x, sensor_data_as_wsws, '--');
            plot(x, decay, ':')
            title(['AS WSWS, L_{inf} = ' num2str(diff_2D_as_wsws)]);
        end

        subplot(4, 1, 3)
        plot(x, sensor_data_2D_x);
        hold on;
        plot(x, sensor_data_as_wsws_fft, '--');
        plot(x, decay, ':')
        title(['AS WSWS-FFT, L_{inf} = ' num2str(diff_2D_as_wsws_fft)]);        

        if exist('dtt1D', 'file')
            subplot(4, 1, 4)
            plot(x, sensor_data_2D_x);
            hold on;
            plot(x, sensor_data_as_wswa, '--');
            plot(x, decay, ':')
            title(['AS WSWA, L_{inf} = ' num2str(diff_2D_as_wswa)]);
        end
        
        scaleFig(1, 1.5);
        drawnow;
        
    end    
    
end