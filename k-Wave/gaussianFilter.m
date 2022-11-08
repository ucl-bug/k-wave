function signal = gaussianFilter(signal, Fs, freq, bandwidth, plot_filter)
%GAUSSIANFILTER Filter signals using a frequency domain Gaussian filter.
%
% DESCRIPTION:
%     gaussianFilter applies a frequency domain Gaussian filter with the
%     specified center frequency and percentage bandwidth to the input
%     signal. If the input signal is given as a matrix, the filter is
%     applied to each matrix row. 
%
% USAGE:
%     signal = gaussianFilter(signal, Fs, freq, bandwidth)
%     signal = gaussianFilter(signal, Fs, freq, bandwidth, plot_filter)
%
% INPUTS:
%     signal      - signal/s to filter
%     Fs          - sampling frequency [Hz]
%     freq        - filter center frequency [Hz]
%     bandwidth   - filter bandwidth [%]
%
% OPTIONAL INPUTS:
%     plot_filter - Boolean controlling whether the filtering process is
%                   plotted
%
% OUTPUTS:
%     signal      - filtered signal/s
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 10th December 2011
%     last update - 4th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby
%
% See also applyFilter, filterTimeSeries, gaussian

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

% suppress m-lint warnings for unused assigned values
%#ok<*ASGLU>

% compute the double-sided frequency axis
N = length(signal(1, :));
if mod(N, 2)==0
    
    % N is even
    f = (-N/2:N/2-1) * Fs/N;
    
else
    
    % N is odd
    f = (-(N-1)/2:(N-1)/2) * Fs/N;
    
end

% compute gaussian filter coefficients
mean = freq;
variance = (bandwidth / 100 * freq / (2 * sqrt(2 * log(2)))).^2;
magnitude = 1;

% create the double-sided Gaussian filter
gauss_filter = max(gaussian(f, magnitude, mean, variance), gaussian(f, magnitude, -mean, variance));

% store a copy of the central signal element if required for plotting
if nargin == 5 && plot_filter
    example_signal = signal(ceil(end/2), :);
end

% apply filter
signal = real(ifft(ifftshift(...
            bsxfun(@times, gauss_filter, ...
                fftshift(fft(signal, [], 2), 2) ...
            )...
        , 2), [], 2));

% plot filter if required
if nargin == 5 && plot_filter
    
    % compute amplitude spectrum of central signal element
    as = fftshift(abs(fft(example_signal)) / N);
    as_filtered = fftshift(abs(fft(signal(ceil(end/2), :))) / N);
    
    % get axis scaling factors
    [f_sc, f_scale, f_prefix] = scaleSI(f);
    
    % produce plot
    figure;
    plot(f .* f_scale, as ./ max(as(:)), 'k-');
    hold on;
    plot(f .* f_scale, gauss_filter, 'b-');
    plot(f .* f_scale, as_filtered ./ max(as(:)), 'r-');
    xlabel(['Frequency [' f_prefix 'Hz]']);
    ylabel('Normalised Amplitude');
    legend('Original Signal', 'Gaussian Filter', 'Filtered Signal');
    axis tight;
    set(gca, 'XLim', [0, f(end) .* f_scale]);
    
end