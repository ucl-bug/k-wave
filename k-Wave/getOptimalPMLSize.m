function pml_sz = getOptimalPMLSize(grid_sz, pml_range, axisymmetric)
%GETOPTIMALPMLSIZE Find PML size to give the smallest prime factors.
%
% DESCRIPTION:
%     getOptimalPMLSize finds the size of the perfectly matched layer (PML)
%     that gives an overall grid size with the smallest prime factors when 
%     using the first-order simulation functions in k-Wave with the
%     optional input 'PMLInside', false. Choosing grid sizes with small
%     prime factors can have a significant impact on the computational
%     speed, as the code computes spatial gradients using the fast Fourier
%     transform (FFT).
%
% USAGE:
%     pml_sz = getOptimalPMLSize(grid_sz)
%     pml_sz = getOptimalPMLSize(grid_sz, pml_range)
%     pml_sz = getOptimalPMLSize(grid_sz, [], axisymmetric)
%     pml_sz = getOptimalPMLSize(grid_sz, pml_range, axisymmetric)
%
% INPUTS:
%     grid_sz         - Grid size defined as a one (1D), two (2D), or three
%                       (3D) element vector. Alternatively, can be an
%                       object of the kWaveGrid class defining the
%                       Cartesian and k-space grid fields.
%
% OPTIONAL INPUTS
%     pml_range       - Two element vector specifying the minimum and
%                       maximum PML size (default = [10, 40]).
%     axisymmetric    - If using the axisymmetric code, string specifying
%                       the radial symmetry. Allowable inputs are 'WSWA'
%                       and 'WSWS' (default = ''). This is important as the
%                       axisymmetric code only applies to the PML to the
%                       outside edge in the radial dimension.
%
% OUTPUTS:
%     pml_opt         - PML size that gives the overall grid with the
%                       smallest prime factors.
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 31st January 2018
%     last update     - 7th January 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby
%
% See also fft

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

% check if grid size is given as kgrid, and extract grid size
if isa(grid_sz, 'kWaveGrid')
    switch grid_sz.dim
        case 1
            grid_sz = grid_sz.Nx;
        case 2
            grid_sz = [grid_sz.Nx, grid_sz.Ny];
        case 3
            grid_sz = [grid_sz.Nx, grid_sz.Ny, grid_sz.Nz];
    end
end

% assign grid size
grid_dim = length(grid_sz);

% check grid size is 1, 2, or 3
if (grid_dim < 1) || (grid_dim > 3)
    error('Grid dimensions must be given as a 1, 2, or 3 element vector.');
end

% check for pml_range input
if (nargin < 2) || (isempty(pml_range))
    pml_range = [10, 40]; 
end

% force integer
pml_range = round(pml_range);

% check for positive values
if any(pml_range < 0)
    error('Optional input pml_range must be positive.');
end

% check for correct length
if length(pml_range) ~= 2
    error('Optional input pml_range must be a two element vector.');
end

% check for monotonic
if pml_range(2) <= pml_range(1)
    error('The second value for pml_range must be greater than the first.');
end

% check for axisymmetric input
if (nargin < 3) || (isempty(axisymmetric))
    axisymmetric = false;
end

% check for correct string
if ischar(axisymmetric) && ~(strcmp(axisymmetric, 'WSWA') || strcmp(axisymmetric, 'WSWS'))
    error('Optional input axisymmetric must be set to ''WSWA'' or ''WSWS''.');
end

% check for correct dimensions
if ischar(axisymmetric) && (grid_dim ~= 2)
    error('Optional input axisymmetric is only valid for 2D grid sizes.');
end

% create array of PML values to search
pml_size = pml_range(1):pml_range(2);

% extract largest prime factor for each dimension for each pml size
facs = zeros(grid_dim, length(pml_size));
for dim = 1:grid_dim
    for index = 1:length(pml_size)
        if ischar(axisymmetric) && (dim == 2)
            switch axisymmetric
                case 'WSWA'
                    facs(dim, index) = max(factor((grid_sz(dim) + pml_size(index)) * 4));
                case 'WSWS'
                    facs(dim, index) = max(factor((grid_sz(dim) + pml_size(index)) * 2 - 2));
            end
        else
            facs(dim, index) = max(factor(grid_sz(dim) + 2 * pml_size(index)));
        end
    end
end

% get best dimension size
[~, ind_opt] = min(facs, [], 2);

% assign output
pml_sz = zeros(1, grid_dim);
pml_sz(1) = pml_size(ind_opt(1));
if grid_dim > 1
    pml_sz(2) = pml_size(ind_opt(2));
end
if grid_dim > 2
    pml_sz(3) = pml_size(ind_opt(3));
end