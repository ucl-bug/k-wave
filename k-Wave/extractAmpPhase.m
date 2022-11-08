function varargout = extractAmpPhase(data, Fs, source_freq, varargin)
%EXTRACTAMPPHASE Extract amplitude and phase from CW signals.
%
% DESCRIPTION:
%     extractAmpPhase extracts the amplitude and phase information at a
%     specified frequency from a vector or matrix of time series data. By
%     default the time dimension is set to the highest non-singleton
%     dimension. The amplitude and phase are extracted from the frequency
%     spectrum, which is calculated using a windowed and zero padded FFT.
%     The values are extracted at the frequency closest to source_freq.
%
% USAGE:
%     amp = extractAmpPhase(data, Fs, source_freq)
%     amp = extractAmpPhase(data, Fs, source_freq, ...)
%     [amp, phase] = extractAmpPhase(data, Fs, source_freq)
%     [amp, phase] = extractAmpPhase(data, Fs, source_freq, ...)
%     [amp, phase, freq] = extractAmpPhase(data, Fs, source_freq)
%     [amp, phase, freq] = extractAmpPhase(data, Fs, source_freq, ...)
%
% INPUTS:
%     data         - matrix of time signals [s]
%     Fs           - sampling frequency [Hz]
%     source_freq  - frequency at which the amplitude and phase should be
%                    extracted [Hz]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Dim'        - dimension over which the signals vary in time
%                    (default = highest non-singleton dimension)
%     'FFTPadding' - scaling parameter used to zero pad the FFT, where
%                    the FFT length = FFTPadding * Nt (default = 3)
%     'Window'     - parameter string controlling the window type used to
%                    filter the signal before the FFT is taken (default =
%                    'Hanning'). Any valid input types for getWin may be
%                    used. 
%
% OUTPUTS:
%     amp          - amplitude [au]
%     phase        - phase between 0 and 2*pi [rad]
%     freq         - closest frequency to source_freq at which amplitude
%                    and phase are extracted [Hz]
%
% ABOUT:
%     author       - Bradley Treeby and Yan To Ling
%     date         - 26th February 2015
%     last update  - 19th February 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2019 Bradley Treeby and Yan To Ling
%
% See also createCWSignals, toneburst

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
dim                     = 'auto';
fft_padding             = 3;
num_req_input_variables = 3;
window                  = 'Hanning';

% check size of input data
if ndims(data) > 4
    error('Input data must have 1, 2, 3, or 4 dimensions.');
end

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Dim'
                dim = varargin{input_index + 1};            
            case 'FFTPadding'
                fft_padding = varargin{input_index + 1};
            case 'Window'
                window = varargin{input_index + 1};
            otherwise
                error('Unknown optional input parameter.');
        end
    end
end

% check for the dim input
if strcmp(dim, 'auto')
    dim = ndims(data);
    if (dim == 2) && (size(data, 2) == 1)
        dim = 1;
    end
end

% create 1D window and reshape to be oriented in the time dimension of the
% input data
[win, coherent_gain] = getWin(size(data, dim), window);
win = reshape(win, [ones(1, dim - 1), length(win), 1]);

% apply window to time dimension of input data
data = bsxfun(@times, win, data);

% compute amplitude and phase spectra
[f, func_as, func_ps] = spect(data, Fs, 'FFTLength', fft_padding .* size(data, dim), 'Dim', dim);

% correct for coherent gain
func_as = func_as ./ coherent_gain;

% find the index of the frequency component closest to source_freq
[~, f_index] = findClosest(f, source_freq);

% get size of output variable, collapsing the time dimension
sz = size(data);
sz(dim) = 1;

% extract amplitude and relative phase at freq_index
amp   = zeros(sz); %#ok<PREALL>
phase = zeros(sz); %#ok<PREALL>

% extract values
switch dim
    case 1
        amp   = func_as(f_index, :, :, :);
        phase = func_ps(f_index, :, :, :);
    case 2
        amp   = func_as(:, f_index, :, :);
        phase = func_ps(:, f_index, :, :);
    case 3
        amp   = func_as(:, :, f_index, :);
        phase = func_ps(:, :, f_index, :);
    case 4
        amp   = func_as(:, :, :, f_index);
        phase = func_ps(:, :, :, f_index);
    otherwise
        error('dim must be 1, 2, 3, or 4');
end

% assign amplitude output
varargout(1) = {amp};

% assign phase output
if nargout > 1
    varargout(2) = {phase};
end

% assign frequency output
if nargout == 3
    varargout(3) = {f(f_index)};    
end