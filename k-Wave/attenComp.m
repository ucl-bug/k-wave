function [signal, tfd, cutoff_freq] = attenComp(signal, dt, c, alpha_0, y, varargin)
%ATTENCOMPEN Attenuation compensation using time variant filtering.
%
% DESCRIPTION:
%     attenComp corrects for frequency dependent acoustic attenuation in
%     photoacoustic signals using time-variant filtering [1]. The time-
%     variant filter is constructed to correct for acoustic attenuation and
%     dispersion following a frequency power law of the form alpha_0*f^y
%     under the assumption the distribution of attenuation parameters is
%     homogeneous. The filter is applied directly to the recorded
%     time-domain signals using a form of non-stationary convolution. The
%     approach is computationally efficient and can be used with any
%     detector geometry or reconstruction algorithm.
%
%     To prevent high-frequency noise from being amplified, the
%     compensation is regularised using a Tukey window with a time-variant
%     cutoff frequency. The cutoff frequency can be specified manually
%     using the optional input 'FilterCutoff'. This is set as a two-element
%     vector corresponding to the cutoff frequency in Hz for the first and
%     last time points, respectively. For a fixed cutoff, these should be
%     specified as the same value, e.g., [3e6, 3e6]. Alternatively, if
%     'FilterCutoff' is set to 'auto' (the default), the cutoff frequency
%     is chosen based on the local time-frequency distribution of the
%     recorded signals using the following steps:          
%
%         1. Compute the average time-frequency distribution of the input
%         signals (the method can be defined using the optional input
%         'Distribution')  
%
%         2. Threshold the time-frequency distribution to remove noise (the
%         threshold value can be defined using the optional input 
%         'NoiseThreshold')   
%
%         3. Calculate the integral of the thresholded time-frequency
%         distribution at each time point using cumsum  
%
%         4. Find the cutoff frequency at each time point where the
%         integral reaches a given percentage of the maximum value (this
%         percentage can be defined using the optional input
%         'EnergyThreshold')   
%
%         5. Increase the filter cutoff frequency by a fixed multiplier so
%         the cutoff corresponds to the edge of the passband for the Tukey
%         window (the multiplier can be defined using the optional input
%         'FrequencyMultiplier')
%
%         6. Smooth the variation of the cutoff frequency over time (the
%         smoothing function can be defined using the optional input
%         'FitType')  
%
%         7. Threshold any values of the cutoff frequency below zero or
%         above the Nyquist limit 
%
%     If the input contains a matrix of signals, the cutoff frequency is
%     based on the average time frequency distribution. To calculate the
%     cutoff frequency for each signal individually, this function should
%     be called in a loop. This can be parallelised, for example, using
%     parfor from the parallel computing toolbox. For further details
%     about this function and attenuation compensation using time variant
%     filtering, see the reference below.      
%
%     [1] B. E. Treeby (2013) "Acoustic attenuation compensation in
%     photoacoustic tomography using time-variant filtering," J. Biomed.
%     Opt., vol. 18, no. 3, p.036008.  
%
% USAGE:
%     signal_comp = attenComp(signal, dt, c, alpha_0, y)
%     signal_comp = attenComp(signal, dt, c, alpha_0, y, ...)
%     [signal_comp, tfd, cutoff_freq] = attenComp(signal, dt, c, alpha_0, y)
%     [signal_comp, tfd, cutoff_freq] = attenComp(signal, dt, c, alpha_0, y, ...)
%
% INPUTS:
%     signal          - time series to compensate, indexed as
%                       (sensor_index, time_index)
%     dt              - time step [s]
%     c               - sound speed [m/s]
%     alpha_0         - power law absorption prefactor [dB/(MHz^y cm)]
%     y               - power law absorption exponent [0 < y < 3, y ~= 1]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings.
%
%     'DisplayUpdates' - Boolean controlling whether command line updates
%                       and compute time are printed to the command line
%                       (default = true).
%     'Distribution'  - Time-frequency distribution used to automatically
%                       compute the filter cutoff frequency if
%                       'FilterCutoff' is set to 'auto'. Valid inputs
%                       are: 
%                           'Rihaczek'  (default)
%                           'Wigner'    
%                       The second option requires the Discrete TFD
%                       toolbox from http://tfd.sourceforge.net
%     'EnergyThreshold' - Threshold value given as a percentage of the
%                       total amplitude spectrum used to choose the
%                       filter cutoff frequency at each time point 
%                       (default = 0.98).
%     'FilterCutoff'  - Option to manually define the cutoff frequencies
%                       for a linear variation in the filter cutoff
%                       instead of using an automatic search. Valid
%                       inputs are: 
%                           'auto'      (default)
%                           [d_min_cutoff, d_max_cutoff]
%     'FitType'       - Fitting type used to smooth the filter cutoff
%                       frequency after an automatic search. Valid inputs
%                       are:  
%                           'spline'    (smoothing spline, the default)
%                           'mav'       (moving average)
%                           'linear'    (linear fit)
%     'FrequencyMultiplier' - By default, the compensation is regularised
%                       using a Tukey window with a time-variant cutoff
%                       frequency. The default Tukey window has a taper
%                       ratio of 0.5, so the filter cutoff frequency
%                       found by the automatic search is increased by a
%                       frequency multiplier so that the filter cutoff
%                       frequency corresponds to the edge of the passband
%                       of the Tukey window (default = 2). 
%     'NumSplines'    - Number of spline segments used in the smoothing
%                       spline if 'FitType' is set to 'spline' 
%                       (default = 40).
%     'NoiseThreshold'- Threshold value given as a percentage of the
%                       signal maximum used to threshold the TFD before
%                       the automatic search for the filter cutoff
%                       (default = 0.03).
%     'Plot'          - Boolean controlling whether a plot of the time
%                       frequency distribution and filter cutoff 
%                       frequency are displayed (default = false). 
%     'PlotRange'     - Option to manually set the plot range in the
%                       frequency axis. Valid options are:
%                           'auto'      (1.5 x maximum filter cutoff)
%                           'full'      (maximum supported frequency range)
%                           [f_min, f_max]
%     'TaperRatio'    - Taper ratio used to construct the Tukey Windows
%                       (default = 0.5).
%     'T0'            - Time index of T0 in the input signals. For
%                       photoacoustic imaging, T0 corresponds to the
%                       arrival of the excitation laser pulse at the
%                       sample (default = 0).
%
% OUTPUTS:
%     signal_comp     - time series after attenuation compensation
%     tfd             - average time frequency distribution of the input
%                       signals
%     cutoff_freq     - filter cutoff frequency for each time index
% 
% ABOUT:
%     author          - Bradley Treeby
%     date            - 15th June 2012
%     last update     - 21st July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2021 Bradley Treeby

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

% =========================================================================
% INPUT CHECKING
% =========================================================================

% define literals
PLOT_DNR        = 40;           % dynamic range used in TFD plot [dB]
PLOT_RANGE_MULT = 1.5;          % multiplier used to scale the auto-plot range

% define defaults for optional input settings
num_req_input_variables = 5;
display_updates = false;        % Boolean option to display command line updates 
distribution    = 'Rihaczek';   % default TF distribution
energy_cutoff   = 0.98;         % cutoff frequency as a [%] of total energy
freq_multiplier = 2;            % used to increase the cutoff_f for a smoother filter 
filter_cutoff   = 'auto';       % automatically compute cutoff based on TFD 
fit_type        = 'spline';     % default fit type used for smooth distribution
noise_cutoff    = 0.03;         % [%] of signal max, used to threshold the signals
num_splines     = 40;           % used for fit_type = 'spline'
plot_tfd        = false;        % plot TFD and cutoff
plot_range      = 'auto';       % plot range
t0              = 1;            % index of laser pulse
taper_ratio     = 0.5;          % taper ratio used for Tukey Window

% NOTES: in JBO paper, shepp logan example uses noise/multp of 0.02/2.5,
% in vivo example uses 0.03/2.

% rotate input signal from (sensor_index, time_index) to (time_index,
% sensor_index)
signal = signal.';

% extract signal characteristics
[N, num_signals] = size(signal);

% convert absorption coefficient to nepers
alpha_0 = db2neper(alpha_0, y);

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Distribution'
                distribution = varargin{input_index + 1};
            case 'DisplayUpdates'
                display_updates = varargin{input_index + 1};
            case 'EnergyThreshold'
                energy_cutoff = varargin{input_index + 1};
            case 'FilterCutoff'
                filter_cutoff = varargin{input_index + 1};
            case 'FitType'             
                fit_type = varargin{input_index + 1};
            case 'FrequencyMultiplier'
                freq_multiplier = varargin{input_index + 1};
            case 'NoiseThreshold'
                noise_cutoff = varargin{input_index + 1};
            case 'NumSplines'
                num_splines = varargin{input_index + 1};
            case 'Plot'
                plot_tfd = varargin{input_index + 1};
            case 'PlotRange'
                plot_range = varargin{input_index + 1};
            case 'TaperRatio'
                taper_ratio = varargin{input_index + 1};
            case 'T0'
                t0 = varargin{input_index + 1};
            otherwise
                error(['Unknown optional input ' varargin{input_index} '.']);
        end
    end
end

% update command line status
if display_updates
    disp('Applying time variant filter...');
end

% check FitType input
if strcmp(fit_type, 'mav')
    
    % define settings for moving average based on the
    % length of the input signals
    mav_terms = round(N*1e-2);
    mav_terms = mav_terms + rem(mav_terms, 2);
    
elseif ~(strcmp(fit_type, 'linear') || strcmp(fit_type, 'spline'))
    
    % throw error for unknown input
    error('Optional input ''FitType'' must be set to ''spline'', ''mav'', or ''linear''.');
    
end

% =========================================================================
% COMPUTE AVERAGE TIME FREQUENCY DISTRIBUTION (TFD) OF INPUT
% =========================================================================

% sample frequency
Fs = 1 / dt;

% create time and frequency axes
t_array = dt * (0:N - 1);

% compute the double-sided frequency axis
if mod(N, 2)==0
    
    % N is even
    f_array = (-N/2:N/2 - 1) * Fs / N;
    
else
    
    % N is odd
    f_array = (-(N - 1)/2:(N - 1)/2) * Fs / N;
    
end

% compute the TFD if required
if strcmp(filter_cutoff, 'auto') || plot_tfd || nargout > 1

    % update display
    if display_updates
        if num_signals > 1
            tic, fprintf('  calculating average time-frequency spectrum...');
        else
            tic, fprintf('  calculating time-frequency spectrum...');
        end
    end

    % compute the TFD of the input signal
    switch distribution
        case 'Rihaczek'

            % loop through signals and calculate average TFD
            tfd = (conj(fft(signal(:, 1))) * signal(:, 1).');
            if num_signals > 1
                for index = 2:num_signals
                    tfd = tfd + (conj(fft(signal(:, index))) * signal(:, index).');
                end
            end
            tfd = tfd .* exp( -1i * (0:2*pi/N:2*pi*(1-1/N))' * (0:N-1));
            tfd = fftshift(tfd, 1) ./ (N*num_signals);

        case 'Wigner'

            % loop through signals and calculate average TFD (this is currently
            % not very optimised and much slower than Rihaczek)
            tfd = qwigner2(signal(:, 1), Fs);
            if num_signals > 1
                for index = 2:num_signals
                    tfd = tfd + qwigner2(signal(:, index), Fs);
                end
            end
            tfd = tfd ./ num_signals;

    end

    % take the absolute value 
    % tfd = abs(real(tfd));
    tfd = abs(tfd);
    
end

% plot the time-frequency spectra of the time signal
if plot_tfd
    
    figure;
    imagesc(t_array .* 1e6, f_array .* 1e-6 , 20 * log10(tfd ./ max(tfd(:))));
    set(gca, 'FontSize', 12, 'YDir', 'normal', 'YLim', [0, f_array(end) * 1e-6]);
    xlabel('Time [\mus]');
    ylabel('Frequency [MHz]');
    if num_signals > 1
        title('Average Time-Frequency Distribution of Signal');
    else
        title('Time-Frequency Distribution of Signal');
    end
    h = colorbar;
    db_plot_scale = [-PLOT_DNR, 0];
    caxis(db_plot_scale);
    xlabel(h, '[dB]');
    
end

% =========================================================================
% FIND CUTOFF FREQUENCIES
% =========================================================================

% calculate filter cutoff frequencies automatically based on TFD
if strcmp(filter_cutoff, 'auto')

    % update display
    if display_updates
        toc, tic, fprintf('  finding filter thresholds... ');
    end

    % shrink t-f_array spectrum to half frequency space to allow cumsum
    f_array_hs = f_array(f_array >= 0);
    tfd_hs = 2 * tfd(f_array >= 0, :);

    % threshold the t-f_array spectrum to remove noise
    threshold = noise_cutoff * max(tfd_hs(:));
    tfd_hs(tfd_hs < threshold) = 0;

    % create integral variables of thresholded t-f_array spectrum
    tfd_int = cumsum(tfd_hs, 1);
    tfd_tot = sum(tfd_hs, 1);

    % preallocate variable
    cutoff_index = zeros(size(t_array));

    % find the cutoff frequency for each time
    for tw_index = 1:length(tfd_tot)  
        [~, cutoff_index(tw_index)] = findClosest(tfd_int(:, tw_index) ./ tfd_tot(tw_index), energy_cutoff);    
    end

    % convert the cutoff index to a frequency and increase by the multiplier
    cutoff_freq_array = (f_array_hs(cutoff_index)) * freq_multiplier;

    % plot the threshold
    if plot_tfd
        hold on;
        plot((0:N - 1) * dt * 1e6, cutoff_freq_array * 1e-6, 'y-');
        plot((0:N - 1) * dt * 1e6, -cutoff_freq_array * 1e-6, 'y-');
        drawnow;
    end

    % smooth the variation in the cutoff frequency from point to point
    switch fit_type

        case 'linear'

            % fit a straight line
            neg_penalty = 10;
            opt_vals = fminsearch(@(x) constlinfit(1:N, cutoff_freq_array, x(1), x(2), neg_penalty), [-f_array_hs(end) / N, f_array_hs(floor(end/2))]);
            cutoff_freq_array = opt_vals(1) * (1:N) + opt_vals(2);

        case 'spline'

            % fit a smoothed spline
            pp = splinefit(1:N, cutoff_freq_array, num_splines, 'r');
            cutoff_freq_array = ppval(pp, 1:N);

        case 'mav'

            % fit a moving average
            cutoff_freq_array = filter(ones(1, mav_terms) * 1/mav_terms, 1, cutoff_freq_array);
            cutoff_freq_array = circshift(cutoff_freq_array, [0, -mav_terms/2]);

        case 'off'

        otherwise
            error('unknown fit_type setting');

    end
else
    
    % set manual values
    cutoff_freq_array = (filter_cutoff(2) - filter_cutoff(1))./(N - 1)*(1:N) + filter_cutoff(1);
    
end

% constrain values outside frequency range
cutoff_freq_array(cutoff_freq_array > f_array(end)) = f_array(end);

% constrain negative values
cutoff_freq_array(cutoff_freq_array < 0) = 0;

% plot final threshold
if plot_tfd
    
    % plot
    hold on;
    plot(t_array * 1e6, cutoff_freq_array * 1e-6, 'w--');
    plot(t_array * 1e6, -cutoff_freq_array * 1e-6, 'w--');

    % set the plot range
    if isnumeric(plot_range)
        set(gca, 'YLim', plot_range);
    elseif strcmp(plot_range, 'auto') && (PLOT_RANGE_MULT * max(cutoff_freq_array(:)) < f_array(end))
        set(gca, 'YLim', PLOT_RANGE_MULT * max(cutoff_freq_array(:)) * [-1, 1] * 1e-6);
    end
    drawnow;
    
end

% =========================================================================
% CREATE FILTER
% =========================================================================

% update display
if display_updates
    toc, tic, fprintf('  creating time variant filter... ');
end

% create distance vector accounting for the index of the laser pulse
% relative to the start of the signals
dist_vec = c * dt * (0 - t0:N - 1 - t0);
dist_vec(dist_vec < 0 ) = 0;

% create the time variant filter
[f_mat, dist_mat] = meshgrid(f_array, dist_vec);
tv_filter = alpha_0 * dist_mat .* ( ...
    (2 * pi * abs(f_mat)).^y - 1i * tan(pi * y / 2) * (2 * pi * f_mat) .* (2 * pi * abs(f_mat)).^(y - 1) ...
    );

% convert cutoff frequency to a window size
N_win_array = floor((cutoff_freq_array / f_array(end)) * N) - 1;
N_win_array(rem(N_win_array, 2) == 1) = N_win_array(rem(N_win_array, 2) == 1) + 1;

% loop through each time/distance and create a row of the filter
for t_index = 1:N
    
    % get the filter cutoff freq
    cutoff_freq = cutoff_freq_array(t_index);
    N_win = N_win_array(t_index);
   
    if cutoff_freq ~= 0

        % create window
        win = zeros(N, 1);
        win_pos = ceil((N - N_win)/2) + 1;
        win(win_pos:win_pos + N_win - 1) = getWin(N_win, 'Tukey', 'Param', taper_ratio, 'Rotation', true);

        % window row of tv_filter
        tv_filter(t_index, :) = tv_filter(t_index, :) .* win.';
    
    else
        tv_filter(t_index, :) = 0;
    end
    
end

% take exponential
tv_filter = exp(tv_filter);

% shift frequency axis (filter_mat is created with 0 in centre) and then
% compute the inverse FFT 
tv_filter = real(ifft(ifftshift(tv_filter, 2), [], 2));

% apply circular shift
for t_index = 1:N
    tv_filter(t_index, :) = circshift(tv_filter(t_index, :), [0, t_index - 1]);
end

% zero out lower and upper triangles
ones_mat = ones(N, N);
tv_filter(tril(ones_mat, -ceil(N/2) + 1) == 1) = 0;
tv_filter(triu(ones_mat, ceil(N/2) + 1) == 1) = 0;

% update display
if display_updates
    toc, fprintf('  applying time variant filter... '), tic;
end

% apply the filter
for index = 1:num_signals
    signal(:, index) = tv_filter * signal(:, index);
end

% rotate output signal from (time_index, sensor_index) back to
% (sensor_index, time_index) 
signal = signal.';

% update display
if display_updates
    toc;
end

% =========================================================================
% FITTING FUNCTION
% =========================================================================

function error = constlinfit(x, y, a, b, neg_penalty)
% function for fminbnd to fit a linear slope with negative values penalised
% Bradley Treeby
% 15th June 2012

% set default value for negative penalty
if nargin < 5
    neg_penalty = 10;
end

% calculate error
error = (a .* x + b) - y;

% penalise negative values
error(error < 0) = error(error < 0) .* neg_penalty;

% give a single error value
error = sum(abs(error(:)));