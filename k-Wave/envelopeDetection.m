function env = envelopeDetection(x)
%ENVELOPEDETECTION Extract signal envelope using the Hilbert Transform.
%
% DESCRIPTION:
%     envelopeDetection applies the Hilbert transform to extract the
%     envelope from an input vector x. If x is a matrix, the envelope along
%     each row is returned. 
%
%     Example:
%         x = toneBurst(10e6, 0.5e6, 10);
%         plot(0:length(x) - 1, x, 'k-', 0:length(x) - 1, envelopeDetection(x), 'r-');
%         legend('Input Signal', 'Envelope');
%
% USAGE:
%     env = envelopeDetection(x)
%
% INPUTS:
%     x           - input function
%
% OUTPUTS:
%     env         - envelope of input function
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 20th December 2010
%     last update - 6th June 2017
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

% if the input is one-dimensional, force to be 1 x Nx
if numDim(x) == 1
    x = reshape(x, 1, []);
end

% compute the FFT of the input function (use zero padding to prevent
% effects of wrapping at the beginning of the envelope), if x is a matrix
% the fft's are computed across each row
X = fft(x, length(x(1, :)) * 2, 2);

% multiply the fft by -1i
X = -1i * X;

% set the DC frequency to zero
X(1) = 0;

% calculate where the negative frequencies start in the FFT
neg_f_index = ceil(length(X(1, :)) / 2) + 1;

% multiply the negative frequency components by -1
X(:, neg_f_index:end) = -1 * X(:, neg_f_index:end);

% compute the Hilbert transform using the inverse fft
z = ifft(X, [], 2);

% extract the envelope
env = abs(x + 1i * z(:, 1:length(x)));