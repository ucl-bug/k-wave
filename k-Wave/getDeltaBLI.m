function f = getDeltaBLI(Nx, dx, x, x0, include_imag)
%GETDELTABLI Exact BLI of an arbitrarily positioned delta function.     
%
% DESCRIPTION:
%     getDeltaBLI calculates the exact BLI of an arbitrarily positioned
%     delta function. For grid dimensions with an evenly-sampled
%     periodicity, a small Nyquist frequency sinusoid is added. This
%     sinusoid is invisible on grid samples and has zero amplitude when the
%     delta function lies on a grid node. It is important when the
%     evaluation points aren't grid nodes, and when the delta function is
%     off-grid. It serves to ensure conjugate symmetry in the BLI's Fourier
%     transform.
%
%     The functions is an extension of http://www.math.udel.edu/~braun/...
%     M428/Matlab/interpolation/triginterp.m 
%
% USAGE:
%     f = getDeltaBLI(Nx, dx, x, x0)
%
% INPUTS:
%     Nx          - Number of grid points in the relevant Cartesian
%                   direction.
%     dx          - Grid point spacing [m].
%     x           - Coordinates at which the BLI is evaluated [m].
%     x0          - Coordinate at which the BLI is centred [m].
%
% OUTPUTS:
%     f           - Value of the BLI at the specified coordinates.
%
% ABOUT:
%     author      - Elliott Wise
%     date        - 1st May 2018
%     last update - 5th December 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Elliott Wise
%
% See also getBLI

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

% ignore imaginary component of even function by default
if nargin < 5
    include_imag = false;
end

% check whether the grid has even or odd samples per period
iseven = (rem(Nx, 2) == 0);

% compute BLI
if iseven

    % compute periodic sinc function
    f = sin(pi .* (x - x0) ./ dx) ./ (Nx * tan(pi * (x - x0) / (Nx * dx)));

    % correct indeterminate points
    f((x - x0) == 0) = 1;

    % add Nyquist sinusoid to ensure conjugate symmetry
    f = f - sin(pi * x0 / dx) / Nx .* sin(pi * x / dx);
    if include_imag
        f = f + 1i * cos(pi * x0 / dx) / Nx .* sin(pi * x / dx);
    end

else

    % compute periodic sinc function
    f = sin(pi * (x - x0) / dx) ./ (Nx * sin(pi * (x - x0) / (Nx * dx)));

    % correct indeterminate points
    f((x - x0) == 0) = 1;

end