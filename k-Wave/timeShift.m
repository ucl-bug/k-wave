function data = timeShift(data, ~, forward)
%TIMESHIFT Shift time series to and from staggered temporal grid.
%
% DESCRIPTION:
%     This function will be deprecated in a future version of k-Wave. Use
%     fourierShift instead.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 5th September 2013
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2017 Bradley Treeby
%
% See also gradientSpect, kspaceFirstOrder1D, kspaceFirstOrder2D,
% kspaceFirstOrder3D, pstdElastic2D, pstdElastic3D

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

% display warning
disp('WARNING: timeShift will be deprecated in a future version of k-Wave. Please changes usage to fourierShift.')

% pass inputs to fourierShift
if nargin < 3 || forward
    data = fourierShift(data, 1/2);
else
    data = fourierShift(data, -1/2);
end