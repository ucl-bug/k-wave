function data = fourierShift(data, shift, dim)
%FOURIERSHIFT Resample data using a Fourier interpolant.
%
% DESCRIPTION:
%     fourierShift resamples the input data along the dimension dim using a
%     regular grid that is offset by the non-dimensional distance shift.
%     The resampling is performed using a Fourier interpolant. 
%
%     This function can be used to shift the acoustic particle velocity
%     recorded by the first-order simulation functions to the regular
%     (non-staggered) temporal grid by setting shift to 1/2.
%
%     Example:
%
%         % define input function
%         dt      = pi/5;
%         t       = 0:dt:(2*pi - dt);
%         y       = cos(t);
% 
%         % shift by 1/2 a grid point spacing
%         y_shift = fourierShift(y, 1/2);
% 
%         % calculate exact function for reference
%         y_ref   = cos(t + dt/2);
% 
%         % plot comparison
%         plot(t, y, 'k-s', t, y_shift, 'b-s', t, y_ref, 'rx');
%         legend('Original Function', 'Shifted Function', 'Reference', ...
%             'Location', 'north');
%
% USAGE:
%       data = fourierShift(data, shift)
%       data = fourierShift(data, shift, dim)
%
% INPUTS:
%       data        - input data
%       shift       - non-dimensional shift, where 0 is no shift, 1/2 is
%                     for a staggered grid, and 1 is a full grid point
%
% OPTIONAL INPUTS:
%       dim         - dimension over which the signals vary in time
%                    (default = highest non-singleton dimension)
%
% OUTPUTS:
%       data        - shifted data
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 22nd September 2015
%       last update - 22nd August 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Bradley Treeby
%
% See also gradientSpect.

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

% check for the dim input
if nargin < 3
    dim = ndims(data);
    if (dim == 2) && (size(data, 2) == 1)
        dim = 1;
    end
end

% define the discretisation of the temporal dimension such that there is
% always a DC component 
N = size(data, dim);
if rem(N, 2) == 0
    
    % grid dimension has an even number of points
    k_vec = (2*pi) .* ((-N/2:N/2-1)/N);
    
else
    
    % grid dimension has an odd number of points
    k_vec = (2*pi) .* ((-(N-1)/2:(N-1)/2)/N);
    
end

% force middle value to be zero in case 1/N is a recurring number and the
% series doesn't give exactly zero 
k_vec(floor(N/2) + 1) = 0;

% put the wavenumber vector in the correct orientation for use with bsxfun
switch dim
    case 1
        k_vec = reshape(k_vec, [], 1);
    case 2
        k_vec = reshape(k_vec, 1, []);
    case 3
        k_vec = reshape(k_vec, 1, 1, []);
    case 4
        k_vec = reshape(k_vec, 1, 1, 1, []);
    otherwise
        error('Input dim must be 1, 2, 3, or 4.');
end

% shift the input using a Fourier interpolant
data = real(ifft(bsxfun(@times, ifftshift( exp(1i .* k_vec .* shift) ), fft(data, [], dim)), [], dim));