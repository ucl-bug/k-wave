function filtered_signal = filterTimeSeries(kgrid, medium, signal, varargin)
%FILTERTIMESERIES Filter signal using the Kaiser windowing method.
%
% DESCRIPTION:
%     filterTimeSeries filters an input time domain signal using a low pass
%     filter applied by applyFilter with a specified cut-off frequency,
%     stop-band attenuation, and transition bandwidth. It uses the Kaiser
%     Windowing method to design the FIR filter, which can be implemented
%     as either a zero phase or linear phase filter. The cutoff frequency
%     is defined by a minimum number of points per wavelength. A smoothing
%     ramp can also be applied to the beginning of the signal to reduce
%     high frequency transients.
%
% USAGE:
%     filtered_signal = filterTimeSeries(kgrid, medium, signal)
%     filtered_signal = filterTimeSeries(kgrid, medium, signal, ...)
%
% INPUTS:
%     kgrid           - k-Wave grid object returned by kWaveGrid
%     medium          - k-Wave medium structure
%     signal          - the time domain signal to filter
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings.
%
%     'PlotSignals'   - Boolean controlling whether the time signal is
%                       displayed before and after filtering (default =
%                       false)
%     'PlotSpectrums' - Boolean controlling whether the amplitude spectrum
%                       is displayed before and after filtering 
%                       (default = false)
%     'PPW'           - the number of points per wavelength used to compute
%                       the filter cutoff frequency, setting to 0 turns off
%                       the filtering (default = 3)  
%     'RampPPW'       - the number of points per wavelength used to compute
%                       the length of the cosine start-up ramp, setting to
%                       0 turns off the start-up ramp (default = 0)
%     'StopBandAtten' - attenuation in decibels in the filter stop band
%                       (default = 60)
%     'TransitionWidth' 
%                     - size of the transition relative to the temporal
%                       sampling frequency (default = 0.1)
%     'ZeroPhase'     - Boolean controlling whether a causal or zero phase
%                       filter is applied (default = false)
%
% OUTPUTS:
%     filtered_signal - the filtered time signal
%
% ABOUT:
%     author          - Bradley Treeby and Ben Cox
%     date            - 3rd December 2009
%     last update     - 25th November 2020
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2020 Bradley Treeby and Ben Cox
%
% See also smooth, spectrum

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

%#ok<*ASGLU>

% number of input variables
num_req_input_variables = 3;

% default filter cut-off frequency
points_per_wavelength = 3;

% default ramp length
ramp_points_per_wavelength = 0;

% default settings for the Kaiser window
stop_band_atten = 60;
transition_width = 0.1;
zero_phase = false; 

% default plot settings
plot_signals = false;
plot_spectrums = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}           
            case 'PlotSignals'
                plot_signals = varargin{input_index + 1};
            case 'PlotSpectrums'
                plot_spectrums = varargin{input_index + 1};   
            case 'PPW'
                points_per_wavelength = varargin{input_index + 1};
            case 'RampPPW'
                ramp_points_per_wavelength = varargin{input_index + 1};
                if islogical(ramp_points_per_wavelength) && ramp_points_per_wavelength
                    ramp_points_per_wavelength = points_per_wavelength;
                end
            case 'StopBandAtten'
                stop_band_atten = varargin{input_index + 1}; 
            case 'TransitionWidth'
                transition_width = varargin{input_index + 1}; 
            case 'ZeroPhase'
                zero_phase = varargin{input_index + 1};                 
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check the input is a row vector
if numDim(signal) == 1
    [m, n] = size(signal);
    if n == 1
        signal = signal.';
        rotate_signal = true;
    else
        rotate_signal = false;
    end
else
    error('Input signal must be a vector.');
end

% update the command line status
disp('Filtering input signal...');

% extract the time step
if strcmp(kgrid.t_array, 'auto')
    error('kgrid.t_array must be explicitly defined.');
end

% compute the sampling frequency
Fs = 1 / kgrid.dt;

% extract the minium sound speed
if isfield(medium, 'sound_speed')
    
    % for the fluid code, use medium.sound_speed
    c0 = min(medium.sound_speed(:));
    
elseif (isfield(medium, 'sound_speed_compression') && isfield(medium, 'sound_speed_shear'))
    
    % for the elastic code, combine the shear and compression sound speeds
    % and remove zeros values
    ss = [medium.sound_speed_compression(:); medium.sound_speed_shear(:)];
    ss(ss == 0) = [];
    c0 = min(ss(:));
    
    % cleanup unused variables
    clear ss;
    
else
    error('The input fields medium.sound_speed or medium.sound_speed_compression and medium.sound_speed_shear must be defined.');
end

% extract the maximum supported frequency (two points per wavelength)
f_max = kgrid.k_max * c0 / (2 * pi);

% calculate the filter cut-off frequency
filter_cutoff_f = 2 * f_max / points_per_wavelength;

% calculate the wavelength of the filter cut-off frequency as a number of
% time steps
filter_wavelength = ((2 * pi / filter_cutoff_f) / kgrid.dt);

% filter the signal if required
if points_per_wavelength ~= 0
    filtered_signal = applyFilter(signal, Fs, filter_cutoff_f, 'LowPass',...
        'ZeroPhase', zero_phase, 'StopBandAtten', stop_band_atten,...
        'TransitionWidth', transition_width, 'Plot', plot_spectrums);
end

% add a start-up ramp if required
if ramp_points_per_wavelength ~= 0
    
    % calculate the length of the ramp in time steps
    ramp_length = round(ramp_points_per_wavelength * filter_wavelength / (2 * points_per_wavelength));

    % create the ramp
    ramp = (-cos( (0:(ramp_length - 1)) * pi / ramp_length ) + 1) / 2;
    
    % apply the ramp
    filtered_signal(1:ramp_length) = filtered_signal(1:ramp_length) .* ramp;
    
end

% restore the original vector orientation if modified
if rotate_signal
    filtered_signal = filtered_signal.';
end

% update the command line status
disp(['  maximum frequency supported by kgrid: ' scaleSI(f_max) 'Hz (2 PPW)']);
if points_per_wavelength ~= 0
    disp(['  filter cutoff frequency: ' scaleSI(filter_cutoff_f) 'Hz (' num2str(points_per_wavelength) ' PPW)']);
end
if ramp_points_per_wavelength ~= 0
    disp(['  ramp frequency: ' scaleSI(2 * pi / (2 * ramp_length * kgrid.dt)) 'Hz (' num2str(ramp_points_per_wavelength) ' PPW)']);
end
disp('  computation complete.');

% plot signals if required
if plot_signals
    figure;
    t_array = (0:length(signal) - 1) * kgrid.dt;
    [~, scale, prefix] = scaleSI(t_array(end));
    subplot(2, 1, 1), plot(t_array * scale, signal, 'k-');
    xlabel(['Time [' prefix 's]']);
    ylabel('Signal Amplitude');
    title('Original Signal');
    subplot(2, 1, 2), plot(t_array * scale, filtered_signal, 'k-');
    xlabel(['Time [' prefix 's]']);
    ylabel('Signal Amplitude');
    title('Filtered Signal');
end