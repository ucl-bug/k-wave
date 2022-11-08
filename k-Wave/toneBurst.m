function signal = toneBurst(sample_freq, signal_freq, num_cycles, varargin)
%TONEBURST Create an enveloped single frequency tone burst.
%
% DESCRIPTION:
%     toneBurst creates an enveloped single frequency tone burst for use in
%     ultrasound simulations. If an array is given for the optional input
%     'SignalOffset', a matrix of tone bursts is created where each row
%     corresponds to a tone burst for each value of the 'SignalOffset'. If
%     a value for the optional input 'SignalLength' is  given, the tone
%     burst/s are zero padded to this length (in samples). 
%
% USAGE:
%     signal = toneBurst(sample_freq, signal_freq, num_cycles)
%     signal = toneBurst(sample_freq, signal_freq, num_cycles, ...)
%
% INPUTS:
%     sample_freq     - sampling frequency [Hz]
%     signal_freq     - frequency of the tone burst signal [Hz]
%     num_cycles      - number of sinusoidal oscillations
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Envelope'      - Envelope used to taper the tone burst. Valid inputs
%                       are:
%
%                           'Gaussian' (the default)
%                           'Rectangular'
%                           [num_ring_up_cycles, num_ring_down_cycles]
%
%                       The last option generates a continuous wave signal
%                       with a cosine taper of the specified length at the
%                       beginning and end.
%
%     'Plot'          - Boolean controlling whether the created tone
%                       burst is plotted.
%     'SignalLength'  - Signal length in number of samples, if longer
%                       than the tone burst length, the signal is
%                       appended with zeros.
%     'SignalOffset'  - Signal offset before the tone burst starts in
%                       number of samples.
%
% OUTPUTS:
%     signal          - created tone burst
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 4th December 2009
%     last update     - 13th January 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby
%
% See also gaussian

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

% set usage defaults
num_req_input_variables = 3;
envelope = 'Gaussian';
signal_length = [];
signal_offset = 0;
plot_signal = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Envelope'
                envelope = varargin{input_index + 1};
            case 'Plot'
                plot_signal = varargin{input_index + 1};  
            case 'SignalOffset'
                signal_offset = varargin{input_index + 1};
                signal_offset = round(signal_offset);       % force integer
            case 'SignalLength'
                signal_length = varargin{input_index + 1};
                signal_length = round(signal_length);       % force integer
            otherwise
                error('Unknown optional input.');
        end
    end
end

% calculate the temporal spacing
dt = 1 / sample_freq;

% create the tone burst
tone_length = num_cycles / (signal_freq);
tone_t = 0:dt:tone_length;
tone_burst = sin(2 * pi * signal_freq * tone_t);
tone_index = round(signal_offset);


% check for ring up and ring down input
if numel(envelope) == 2
    
    % assign the inputs
    num_ring_up_cycles = envelope(1);
    num_ring_down_cycles = envelope(2);
    
    % check signal is long enough for ring up and down
    if num_cycles < (num_ring_up_cycles + num_ring_down_cycles)
        error('Input num_cycles must be longer than num_ring_up_cycles + num_ring_down_cycles.');
    end
    
    % get period
    period = 1 ./ signal_freq;
    
    % create x-axis for ramp between 0 and pi
    up_ramp_length_points   = round(num_ring_up_cycles .* period ./ dt);
    down_ramp_length_points = round(num_ring_down_cycles .* period ./ dt);
    up_ramp_axis            = 0:(pi / (up_ramp_length_points - 1)):pi;
    down_ramp_axis          = 0:(pi / (down_ramp_length_points - 1)):pi;
     
    % create ramp using a shifted cosine
    up_ramp                 = (-cos(up_ramp_axis)   + 1) .* 0.5;
    down_ramp               = (cos(down_ramp_axis) + 1) .* 0.5;

    % apply the ramps
    tone_burst(1:up_ramp_length_points) = tone_burst(1:up_ramp_length_points) .* up_ramp;
    tone_burst(end - down_ramp_length_points + 1:end) = tone_burst(end - down_ramp_length_points + 1:end) .* down_ramp;
    
else

    % create the envelope
    switch envelope
        case 'Gaussian'
            x_lim = 3;
            window_x = -x_lim:( 2 * x_lim / (length(tone_burst) - 1) ):x_lim;
            window = gaussian(window_x, 1, 0, 1);
        case 'Rectangular'
            window = ones(size(tone_burst));
        case 'RingUpDown'


        otherwise
            error(['Unknown envelope ' envelope '.']);
    end  

    % apply the envelope
    tone_burst = tone_burst .* window;
    
    % force the ends to be zero by applying a second window
    if strcmp(envelope, 'Gaussian')
        tone_burst = tone_burst .* (getWin(length(tone_burst), 'Tukey', 'Param', 0.05).');
    end
    
end

% % calculate the expected FWHM in the frequency domain
% t_var = tone_length/(2*x_lim);
% w_var = 1/(4*pi^2*t_var);
% fw = 2 * sqrt(2 * log(2) * w_var)

% create the signal with the offset tone burst
if isempty(signal_length)
    signal = zeros(length(tone_index), max(signal_offset(:)) + length(tone_burst));
else
    signal = zeros(length(tone_index), signal_length);
end
for offset = 1:length(tone_index)
    signal(offset, tone_index(offset) + 1:tone_index(offset) + length(tone_burst)) = tone_burst;
end

% plot the signal if required
if plot_signal
    
    % compute suitable axis scaling
    time_axis = (0:length(signal) - 1) * dt;
    [~, scale, prefix] = scaleSI(max(time_axis(:))); 
    
    % create figure
    figure;
    if numel(signal_offset) == 1
        plot(time_axis * scale, signal, 'k-');
    else
        plot(time_axis * scale, signal + repmat(2 * max(signal(:)) * ...
            (1:length(signal_offset)).', 1, length(time_axis)), 'k-');
    end
    xlabel(['Time [' prefix 's]']);
    ylabel('Signal Amplitude');
    axis tight;

end