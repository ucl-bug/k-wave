function filtered_signal = applyFilter(signal, Fs, cutoff_f, filter_type, varargin)
%APPLYFILTER Filter input with low, high, or band pass filter.
%
% DESCRIPTION:
%     applyFilter filters an input signal using filter. The FIR filter
%     coefficients are based on a Kaiser window with the specified cut-off
%     frequency and filter type ('HighPass', 'LowPass' or 'BandPass'). Both
%     causal and zero phase filters can be applied.
%       
% USAGE:
%     filtered_signal = applyFilter(signal, Fs, cutoff_f, filter_type)
%     filtered_signal = applyFilter(signal, Fs, cutoff_f, filter_type, ...)
% 
% INPUTS:
%     signal            - signal to filter
%     Fs                - sampling frequency [Hz]
%     cutoff_f          - filter cutoff frequency/s [Hz]
%     filter_type       - 'HighPass', 'LowPass' or 'BandPass'
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings.
%
%     'Plot'            - Boolean controlling whether the amplitude
%                         spectrum is displayed before and after filtering
%                         (default = false)
%     'StopBandAtten'   - attenuation in decibels in the stop band
%                         (default = 60)
%     'TransitionWidth' - size of the transition based on the temporal
%                         sampling frequency (default = 0.1)
%     'ZeroPhase'       - Boolean controlling whether a zero phase filter
%                         is used (default = false)
%     
% OUTPUTS:
%     filtered_signal   - the filtered signal
%
% ABOUT:
%     author            - Ben Cox and Bradley Treeby
%     date              - 31st December 2009 (previously lowPassFilter)
%     last update       - 21st February 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Ben Cox and Bradley Treeby
%
% See also filter, filterTimeSeries

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

% set optional input defaults
num_req_input_variables     = 4;
zero_phase                  = false;
transition_width            = 0.1;  % as proportion of sampling frequency
stop_band_atten             = 60;   % [dB]
plot_filter                 = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Plot'
                plot_filter = varargin{input_index + 1};
            case 'StopBandAtten'
                stop_band_atten = varargin{input_index + 1};
            case 'TransitionWidth'
                transition_width = varargin{input_index + 1};
            case 'Window'
                w = varargin{input_index + 1};
            case 'ZeroPhase'
                zero_phase = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% store the specified cutoff frequency before modificiation for plotting
if plot_filter
    cutoff_f_plot = cutoff_f;
end

% for a bandpass filter, use applyFilter recursively
if strcmp(filter_type, 'BandPass')
    
    % apply the low pass filter
    func_filt_lp = applyFilter(signal, Fs, cutoff_f(2), 'LowPass', 'StopBandAtten', stop_band_atten, 'TransitionWidth', transition_width, 'ZeroPhase', zero_phase);
    
    % apply the high pass filter
    filtered_signal = applyFilter(func_filt_lp, Fs, cutoff_f(1), 'HighPass', 'StopBandAtten', stop_band_atten, 'TransitionWidth', transition_width, 'ZeroPhase', zero_phase);
    
else
    
    % check filter type
    switch filter_type
        case 'LowPass'
            high_pass = false;
        case 'HighPass'
            high_pass = true;
            cutoff_f = (Fs/2 - cutoff_f);
        otherwise
            error(['Unknown filter type ' filter_type]);
    end

    % make sure input is the correct way around
    [m, n] = size(signal);
    if m > n
        signal = signal.';
    end

    % correct the stopband attenuation if a zero phase filter is being used
    if zero_phase
        stop_band_atten = stop_band_atten/2;
    end

    % decide the filter order
    N = ceil((stop_band_atten - 7.95) / (2.285*(transition_width*pi)));

    % construct impulse response of ideal bandpass filter h(n), a sinc function
    fc = cutoff_f/Fs; % normalised cut-off
    n = -N/2:N/2-1;
    h = 2*fc*sinc(2*pi*fc*n);

    % if no window is given, use a Kaiser window
    if ~exist('w', 'var')

        % compute Kaiser window parameter beta
        if stop_band_atten > 50
            beta = 0.1102*(stop_band_atten - 8.7);
        elseif stop_band_atten >= 21
            beta = 0.5842*(stop_band_atten - 21)^0.4 + 0.07886*(stop_band_atten - 21);
        else
            beta = 0;
        end

        % construct the Kaiser smoothing window w(n)
        m = 0:N-1;
        w = real(besseli(0,pi*beta*sqrt(1-(2*m/N-1).^2)))/real(besseli(0,pi*beta));

    end

    % window the ideal impulse response with Kaiser window to obtain the FIR
    % filter coefficients hw(n) 
    hw = w.*h;

    % modify to make a high_pass filter
    if high_pass
        hw = (-1*ones(1, length(hw))).^(1:length(hw)).*hw;
    end

    % add some zeros to allow the reverse (zero phase) filtering room to work
    L = length(signal);    % length of original input signal
    filtered_signal = [zeros(1,N) signal];

    % apply the filter
    filtered_signal = filter(hw,1,filtered_signal);
    if zero_phase
        filtered_signal = fliplr(filter(hw,1,filtered_signal(L+N:-1:1)));
    end

    % remove the part of the signal corresponding to the added zeros
    filtered_signal = filtered_signal(N+1:L+N);
    
end

% plot the amplitude spectrum and the filter if required
if plot_filter
    
    % calculate the amplitude spectrum of the input signal
    [f, func_as] = spect(signal, Fs);
        
    % plot the input signal
    figure;
    [f_sc, scale, prefix] = scaleSI(max(f));  %#ok<ASGLU>
    plot(f*scale, func_as, 'k-');
    xlabel(['Frequency [' prefix 'Hz]']);
    ylabel('Signal Amplitude');    
    hold on;
    
    % get the axis limits
    xlim = get(gca, 'XLim');
    
    % compute the amplitude spectrum of the filtered signal
    [f, func_as] = spect(filtered_signal, Fs);
    
    % plot the filtered signal
    plot(f*scale, func_as, 'b-');
    set(gca, 'XLim', xlim);
    ylim = get(gca, 'YLim');
   
    % plot the filter cutoff frequency
    if strcmp(filter_type, 'BandPass')
        line([cutoff_f_plot(1)*scale, cutoff_f_plot(1)*scale], [0, ylim(2)], 'LineStyle','--', 'Color', 'k');
        line([cutoff_f_plot(2)*scale, cutoff_f_plot(2)*scale], [0, ylim(2)], 'LineStyle','--', 'Color', 'k');
        legend('input signal', 'filtered signal', 'filter cutoff', 'Location', 'best');         
    else
        line([cutoff_f_plot*scale, cutoff_f_plot*scale], [0, ylim(2)], 'LineStyle','--', 'Color', 'k');
        legend('input signal', 'filtered signal', 'filter cutoff', 'Location', 'best');  
    end

end