function cw_signal = createCWSignals(t_array, freq, amp, phase, ramp_length)
%CREATECWSIGNALS Generate array of CW signals from amplitude and phase.
%
% DESCRIPTION:
%     createCWSignal generates a series of continuous wave (CW) signals
%     based on the 1D or 2D input matrices amp and phase, where each signal
%     is given by:
%
%         amp(i, j) .* sin(2 .* pi .* freq .* t_array + phase(i, j)); 
%
%     To avoid startup transients, a cosine tapered up-ramp is applied to
%     the beginning of the signal. By default, the length of this ramp is 
%     four periods of the wave. The up-ramp can be turned off by setting
%     the ramp_length to 0.
%
%     Example:
%
%         % define sampling parameters
%         f = 5e6;
%         T = 1/f;
%         Fs = 100e6;
%         dt = 1/Fs;
%         t_array = 0:dt:10*T;
%
%         % define amplitude and phase
%         amp = getWin(9, 'Gaussian');
%         phase = linspace(0, 2*pi, 9).';
%
%         % create signals and plot
%         cw_signal = createCWSignals(t_array, f, amp, phase);
%         stackedPlot(cw_signal);
%
% USAGE:
%     cw_signal = createCWSignal(t_array, freq, amp, phase)
%     cw_signal = createCWSignal(t_array, freq, amp, phase, ramp_length)
%
% INPUTS:
%     t_array     - 1D vector of time points [s]
%     freq        - frequency of the CW signal [Hz]
%     amp         - 1D or 2D matrix of amplitudes [au]
%     phase       - 1D or 2D matrix of phases [rad]
%
% OPTIONAL INPUTS:
%     ramp_length - length of the up-ramp used to reduce start-up
%                   transients in periods (default = 4) 
%
% OUTPUTS:
%     cw_signal   - matrix of CW signals
%
% ABOUT:
%     author      - Bradley Treeby and Yan To Ling
%     date        - 4th March 2015
%     last update - 28th January 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2018 Bradley Treeby and Yan To Ling
%
% See also extractAmpPhase, toneburst

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

% expand the phase value if given as a scalar
if numel(phase) == 1
    phase = phase .* ones(size(amp));
end

% check input dimensions
if any(size(amp) ~= size(phase))
    error('Inputs amp and phase must be of equal size.'); 
end

% check for ramp input
if nargin == 4
    ramp_length = 4; 
elseif nargin ~= 5
    error('Incorrect number of inputs.');
end

% get size of input
[N1, N2] = size(amp);

% create input signals
cw_signal = zeros(N1, N2, length(t_array));

% create signal
for index1 = 1:N1
    for index2 = 1:N2
        cw_signal(index1, index2, :) = amp(index1, index2) .* sin(2 .* pi .* freq .* t_array + phase(index1, index2));       
    end
end  

% apply ramp to avoid startup transients
if ramp_length ~= 0

    % get period and time step (assuming dt is constant)
    period = 1 ./ freq;
    dt = t_array(2) - t_array(1);
    
    % create x-axis for ramp between 0 and pi
    ramp_length_points = round(ramp_length .* period ./ dt);
    ramp_axis = 0:(pi / (ramp_length_points - 1)):pi;
     
    % create ramp using a shifted cosine
    ramp = (-cos(ramp_axis) + 1) .* 0.5;
    ramp = reshape(ramp, 1, 1, []);

    % apply ramp to all signals simultaneously
    cw_signal(:, :, 1:ramp_length_points) = bsxfun(@times, ramp, cw_signal(:, :, 1:ramp_length_points));
    
end

% remove singleton dimensions
cw_signal = squeeze(cw_signal);

% if only a single amplitude and phase is given, force time to be the
% second dimensions
if numel(amp) == 1
    cw_signal = reshape(cw_signal, 1, []);
end