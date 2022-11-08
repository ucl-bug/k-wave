function [signal, snr] = addNoise(signal, snr, mode)
%ADDNOISE Add Gaussian noise to a signal for a given SNR.
%
% DESCRIPTION:
%     addNoise adds Gaussian random noise to a one dimensional input signal
%     given the desired signal snr (signal to noise ratio) in decibels. By
%     default, the magnitude of the added noise is calculated based on the
%     rms level of the input signal. For impulsive signals, the optional
%     input mode should be set to 'peak' so the magnitude of the added
%     noise is calculated based on the peak level of the input signal.
%
% USAGE:
%     signal = addNoise(signal, snr)
%     signal = addNoise(signal, snr, mode)
%     [signal, snr] = addNoise(signal, snr)
%     [signal, snr] = addNoise(signal, snr, mode)
%
% INPUTS:
%     signal      - input signal
%     snr         - desired signal snr (signal to noise ratio)
%                   in decibels after adding noise 
%
% OPTIONAL INPUTS:
%     mode        - 'rms' (default) or 'peak'
%
% OUTPUTS:
%     signal      - signal with added noise
%     snr         - actual snr of output signal
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 9th April 2010
%     last update - 21st February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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

if nargin < 3 || strcmp(mode, 'rms')
    
    % calculate the rms amplitude of the input signal
    signal_rms = sqrt(mean(signal(:).^2));
    
elseif strcmp(mode, 'peak')
    
    % calculate the peak amplitude of the input signal
    signal_rms = max(signal(:));
    
else
    
    % throw error if mode input not recognised
    error(['Unknown parameter ' mode ' for optional input mode.']);
    
end

% calculate the standard deviation of the Gaussian noise
std_dev = signal_rms ./ (10.^(snr ./ 20)); 

% calculate noise
noise = std_dev .* randn(size(signal));

% check the snr
noise_rms = sqrt(mean(noise(:).^2));
snr = 20 .* log10(signal_rms ./ noise_rms);

% add noise to the recorded sensor data
signal = signal + noise;