function y = sinc(x)
%SINC Mathematical sinc function.     
%
% DESCRIPTION:
%     sinc calculates the sinc of the input x, where the sinc function is
%     defined as sin(x) / x. This differs from the sinc function in the
%     MATLAB signal processing toolbox which uses sinc(pi*x) / (pi*x).
%
% USAGE:
%     y = sinc(x)
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 14th January 2009
%     last update     - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

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

% warning('off', 'MATLAB:divideByZero');
% y = sin(x)./x;
% y(x == 0) = 1;

% new implementation to avoid changing warning status
zero_vals = (x == 0);
y = sin(x + pi * zero_vals) ./ (x + zero_vals) + zero_vals;