function r = makePixelMap(varargin)
%MAKEPIXELMAP Create matrix of grid point distances from the centre point.
%
% DESCRIPTION:
%     makePixelMap generates a matrix populated with values of how far each
%     pixel in a grid is from the centre (given in pixel coordinates). Both
%     single and double pixel centres can be used by setting the optional
%     input parameter 'OriginSize'. For grids where the dimension size and
%     centre pixel size are not both odd or even, the optional input
%     parameter 'Shift' can be used to control the location of the
%     centerpoint. 
%
%     Examples for a 2D pixel map:
%
%     Single pixel origin size for odd and even (with 'Shift' = [1 1] and 
%     [0 0], respectively) grid sizes:       
%                   
%     x x x       x x x x         x x x x
%     x 0 x       x x x x         x 0 x x
%     x x x       x x 0 x         x x x x
%                 x x x x         x x x x
%
%     Double pixel origin size for even and odd (with 'Shift' = [1 1] and 
%     [0 0], respectively) grid sizes:        
%                   
%     x x x x      x x x x x        x x x x x
%     x 0 0 x      x x x x x        x 0 0 x x
%     x 0 0 x      x x 0 0 x        x 0 0 x x
%     x x x x      x x 0 0 x        x x x x x
%                  x x x x x        x x x x x
%
%     By default a single pixel centre is used which is shifted towards
%     the final row and column.
%
% USAGE:
%     r = makePixelMap(Nx, Ny)
%     r = makePixelMap(Nx, Ny, varargin)
%     r = makePixelMap(Nx, Ny, Nz)
%     r = makePixelMap(Nx, Ny, Nz, varargin)
%
% INPUTS:
%     Nx, Ny, Nz   - size of Cartesian grid in pixels
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings.
%
%     'Shift'      - 2 element binary vector controlling the location of the
%                    centre pixel in the x (row) and y (column) directions
%                    for even grid sizes. A value of 0 will shift towards
%                    the first pixel, and 1 will shift towards the final
%                    pixel (default = [1 1]).
%     'OriginSize' - parameter controlling whether a 'single' or 'double'
%                    pixel is used as the centre (zero) position (default =
%                    'single').  
%
% OUTPUTS:
%     r            - pixel map
%
% ABOUT:
%     author       - Bradley Treeby
%     date         - 6th May 2009
%     last update  - 9th April 2017
% 
% See also makePixelMapPlane, makePixelMapPoint
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

% UPDATES:
%       2009/06/19  - completely overhauled computation and added usage of
%                     a double centre point
%       2009/07/01  - added three-dimensional input options
%       2011/07/19  - re-ordered inputs
%       2017/04/089 - updated documentation and applied coding standards

% define defaults
origin_size = 'single';
shift_def = 1;

% detect whether the inputs are for two or three dimensions
if nargin < 3 || ~isnumeric(varargin{3})
    map_dimension = 2;
    Nx = varargin{1};
    Ny = varargin{2};
    shift = [shift_def, shift_def];
else
    map_dimension = 3;
    Nx = varargin{1};
    Ny = varargin{2};
    Nz = varargin{3};
    shift = [shift_def, shift_def, shift_def];
end

% replace with user defined values if provided
if nargin > map_dimension
    if rem(nargin - map_dimension, 2)
        error('Optional inputs must be entered as param, value pairs.');
    end
    for input_index = map_dimension + 1:2:length(varargin)
        switch varargin{input_index}
            case 'Shift'
                shift = varargin{input_index + 1};
            case 'OriginSize'
                origin_size = varargin{input_index + 1};                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% catch input errors
if ~(strcmp(origin_size, 'single') || strcmp(origin_size, 'double'))
    error('Unknown setting for optional input Center.');
end
if numel(shift) ~= map_dimension
    error(['Optional input Shift must have ' num2str(map_dimension) ' elements for ' num2str(map_dimension) ' dimensional input parameters.']);
end

switch map_dimension
    case 2
        
        % create the maps for each dimension
        nx = createPixelDim(Nx, origin_size, shift(1));
        ny = createPixelDim(Ny, origin_size, shift(2));

        % create plaid grids
        [r_x, r_y] = ndgrid(nx, ny); 

        % extract the pixel radius
        r = sqrt(r_x.^2 + r_y.^2);
        
    case 3
        
        % create the maps for each dimension
        nx = createPixelDim(Nx, origin_size, shift(1));
        ny = createPixelDim(Ny, origin_size, shift(2));
        nz = createPixelDim(Nz, origin_size, shift(3));

        % create plaid grids
        [r_x, r_y, r_z] = ndgrid(nx, ny, nz); 

        % extract the pixel radius
        r = sqrt(r_x.^2 + r_y.^2 + r_z.^2);
        
end

function nx = createPixelDim(Nx, origin_size, shift)
% Nested function to create the pixel radius variable

% grid dimension has an even number of points
if rem(Nx, 2) == 0
    
    % pixel numbering has a single centre point
    if strcmp(origin_size, 'single')
        
        % centre point is shifted towards the final pixel
        if shift == 1
            nx = (-Nx/2:1:Nx/2-1).';
            
        % centre point is shifted towards the first pixel
        else
            nx = (-Nx/2+1:1:Nx/2).';
        end
        
    % pixel numbering has a double centre point    
    else
        nx = [-Nx/2+1:1:0, 0:1:Nx/2-1].';
    end
    
% grid dimension has an odd number of points    
else
    
    % pixel numbering has a single centre point
    if strcmp(origin_size, 'single')
        nx = (-(Nx-1)/2:1:(Nx-1)/2).';
        
    % pixel numbering has a double centre point    
    else

        % centre point is shifted towards the final pixel        
        if shift == 1
            nx = [-(Nx-1)/2:1:0, 0:1:(Nx-1)/2-1].';
            
        % centre point is shifted towards the first pixel            
        else
            nx = [-(Nx-1)/2+1:1:0, 0:1:(Nx-1)/2].';
        end
    end
end