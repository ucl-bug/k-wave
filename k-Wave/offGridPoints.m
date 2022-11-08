function mask = offGridPoints(kgrid, points, scale, varargin)
%OFFGRIDPOINTS Create a non-binary source mask from Cartesian points.
%
% DESCRIPTION:
%     offGridPoints creates a non-binary source mask for defining source or
%     sensor points not aligned with the Cartesian grid points defined by
%     kgrid. For each point given in the input points, a band-limited
%     interpolant (BLI) is computed corresponding to a point source at that
%     location. The point sources are then summed to give the source mask. 
%
%     The input scale can be used to scale the overall mask (e.g., to
%     account for the relative density of the off-grid points compared to
%     the density of kgrid grid points). If given as a vector the same
%     length as points, it can also be used to scale the magnitude of the
%     individual BLIs. If undefined, no scaling is used.
%
%     By default, the spatial extent of the BLI at each point is truncated
%     at 10% of its maximum value. This tolerance can be controlled by
%     using the optional input 'BLITolerance'. If set to 0, the BLI is
%     evaluated at all points on the grid. For low BLI tolerance values, it
%     may be more efficient to set the value to 0, and then truncate the
%     resulting mask.
%
%     The BLI is computed using a sinc function by default, but this is an
%     approximation to the true BLI for a point source. The exact BLI can
%     be used by setting the optional input 'BLIType' to 'exact' (this will
%     also force 'BLITolerance' to be zero).
%
% USAGE:
%     mask = offGridPoints(kgrid, points, scale)
%     mask = offGridPoints(kgrid, points, scale, ...)
%
% INPUTS:
%     kgrid             - Object of the kWaveGrid class defining the
%                         Cartesian and k-space grid fields.
%     points            - List of Cartesian points defined by a matrix with
%                         dimensions num_dims x num_points.
%     scale             - Scaling factor accounting for density of source
%                         points relative to the density of kgrid nodes.
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'BLITolerance'    - Scalar value controlling where the spatial extent
%                         of the BLI at each point is truncated as a
%                         portion of the maximum value (default = 0.1).
%     'BLIType'         - String controlling the BLI expression that is
%                         used for each point source, either 'sinc' or
%                         'exact' (default = 'sinc'). BLITolerance is
%                         ignored if 'exact' is specified.
%     'MaskOnly'        - Boolean controlling whether a logical mask is
%                         returned instead of the non-binary source mask,
%                         where the mask contains the extent of the
%                         off-grid source (defaul = false).
%     'SinglePrecision' - Boolean controlling whether the mask is returned
%                         in single precision. If 'BLITolerance' > 0, then
%                         calculations are also performed in single
%                         precision to improve performance (default =
%                         false).
%     'WaitBar'         - Boolean controlling whether a waitbar is
%                         displayed (default = false).
%
% OUTPUTS:
%     mask              - Non-binary source mask.
%
% ABOUT:
%     author            - Elliott Wise and Bradley Treeby
%     date              - 5th December 2016
%     last update       - 14th July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2016-2021 Elliott Wise and Bradley Treeby
%
% See also kWaveArray

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

% Notes:
% - Scaling effect resulting from truncating the BLI is not accounted for.

% define defaults
bli_tolerance           = 0.1;
bli_type                = 'sinc';
display_wait_bar        = false;
wait_bar_update_freq    = 100;
mask_only               = false;
single_precision        = false;
debug                   = false;

% check dimensions of points input
if size(points, 1) ~= kgrid.dim
    error('Input points must be given as matrix with dimensions num_dims x num_points.');
end

% get the number of off-grid points
num_points = size(points, 2);

% check for scale input
if nargin < 3 || isempty(scale)
    scale = 1;
end

% expand scale value if scalar
if numel(scale) == 1
    scale = scale .* ones(num_points, 1);
elseif numel(scale) ~= num_points
    error('Input scale must be scalar or the same length as points.');
end

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'BLITolerance'
                
                % assign input
                bli_tolerance = varargin{input_index + 1};
                
                % check parameter
                validateattributes(bli_tolerance, {'numeric'}, {'real', 'scalar', '>', 0, '<', 1}, 'offGridPoints', '''BLITolerance''');
                
            case 'BLIType'
                
                % assign input
                bli_type = varargin{input_index + 1};
                
                % check parameter
                validateattributes(bli_type, {'char'}, {}, 'offGridPoints', '''BLIType''');
                if ~any(strcmp(bli_type, {'sinc', 'exact'}))
                    error('Optional input ''BLIType'' must be either ''sinc'' or ''exact''.')
                end
                
            case 'Debug'
                
                % assign input 
                debug = varargin{input_index + 1};
                
                % check parameter
                if ~islogical(debug)
                    error('Optional input ''Debug'' must be Boolean.');
                end                  
                
            case 'MaskOnly'
                
                % assign input
                mask_only = varargin{input_index + 1};
                
                % check parameter
                if ~islogical(mask_only)
                    error('Optional input ''MaskOnly'' must be Boolean.');
                end  
                
            case 'SinglePrecision'
                
                % assign input
                single_precision = varargin{input_index + 1};
                
                % check parameter
                if ~islogical(single_precision)
                    error('Optional input ''SinglePrecision'' must be Boolean.');
                end  
                
            case 'WaitBar'
                
                % assign input
                display_wait_bar = varargin{input_index + 1};
                
                % check parameter
                if ~islogical(display_wait_bar)
                    error('Optional input ''WaitBar'' must be Boolean.');
                end              
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% reassign bli_tolerance if bli_type is 'exact'
if strcmp(bli_type, 'exact')
    bli_tolerance = 0;
end

% extract grid vectors (to avoid recomputing them below)
x_vec = kgrid.x_vec;
y_vec = kgrid.y_vec; 
z_vec = kgrid.z_vec;

% preallocate some variables for speed
if bli_tolerance == 0
    pi_on_dx = pi ./ kgrid.dx;
    pi_on_dy = pi ./ kgrid.dy;
    pi_on_dz = pi ./ kgrid.dz;
else
    scalar_dxyz = true;
    pi_on_dxyz = pi ./ kgrid.dx;
    switch kgrid.dim
        case 2
            if ~(kgrid.dx == kgrid.dy)
                scalar_dxyz = false;
                pi_on_dxyz = pi ./ [kgrid.dx, kgrid.dy];
            end
        case 3
            if ~((kgrid.dx == kgrid.dy) && (kgrid.dx == kgrid.dz))
                scalar_dxyz = false;
                pi_on_dxyz = pi ./ [kgrid.dx, kgrid.dy, kgrid.dz];
            end
    end
end

% initialise the source mask
if mask_only
    mask_type = 'logical';
elseif single_precision
    mask_type = 'single';
else
    mask_type = 'double';
end
switch kgrid.dim
    case 1
        mask = zeros(kgrid.Nx, 1, mask_type);
    case 2
        mask = zeros(kgrid.Nx, kgrid.Ny, mask_type);
    case 3
        mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz, mask_type);
end

% display wait bar
if display_wait_bar
    h = waitbar(0, 'Computing off-grid source mask...');
end

% add to the overall mask using contributions from each source point
for point_ind = 1:num_points
    
    % extract a single point
    point = points(:, point_ind);
    
    % convert to the computational coordinate if the physical coordinate is
    % sampled nonuniformly
    if kgrid.nonuniform
        [point, BLIscale] = mapPoint(kgrid, point);
    end
    
    if bli_tolerance == 0
        
        if mask_only
            mask = true(size(mask));
            return
        else
            switch kgrid.dim
                case 1

                    % evaluate a BLI centered on the point at all the grid nodes
                    % add contribution to overall source mask
                    switch bli_type
                        case 'sinc'
                            mask = mask + scale(point_ind) .* sinc(pi_on_dx * (x_vec - point(1)) );
                        case 'exact'
                            mask = mask + scale(point_ind) .* getDeltaBLI(kgrid.Nx, kgrid.dx, x_vec, point(1));
                    end

                case 2

                    % evaluate a BLI centered on the point at all the grid nodes
                    switch bli_type
                        case 'sinc'
                            mask_t_x = sinc(pi_on_dx * (x_vec - point(1)) );
                            mask_t_y = sinc(pi_on_dy * (y_vec - point(2)) );
                        case 'exact'
                            mask_t_x = getDeltaBLI(kgrid.Nx, kgrid.dx, x_vec, point(1));
                            mask_t_y = getDeltaBLI(kgrid.Ny, kgrid.dy, y_vec, point(2));
                    end

                    % add this contribution to the overall source mask
                    mask = mask + scale(point_ind) .* (mask_t_x * mask_t_y.');

                case 3

                    % evaluate a BLI centered on the point at all the grid nodes
                    switch bli_type
                        case 'sinc'
                            mask_t_x = sinc(pi_on_dx * (x_vec - point(1)) );
                            mask_t_y = sinc(pi_on_dy * (y_vec - point(2)) );
                            mask_t_z = sinc(pi_on_dz * (z_vec - point(3)) );
                        case 'exact'
                            mask_t_x = getDeltaBLI(kgrid.Nx, kgrid.dx, x_vec, point(1));
                            mask_t_y = getDeltaBLI(kgrid.Ny, kgrid.dy, y_vec, point(2));
                            mask_t_z = getDeltaBLI(kgrid.Nz, kgrid.dz, z_vec, point(3));
                    end

                    % add this contribution to the overall source mask
                    mask = mask + scale(point_ind) .* ...
                        reshape(kron(mask_t_y * mask_t_z.', mask_t_x), [kgrid.Nx, kgrid.Ny, kgrid.Nz]);

            end
        end
        
    else
    
        % create an array of neighbouring grid points for BLI evaluation
        switch kgrid.dim
            case 1
                [ind, is] = tolStar(bli_tolerance, kgrid, point, debug);
                xs = x_vec(is);
                xyz = xs;
            case 2
                [ind, is, js] = tolStar(bli_tolerance, kgrid, point, debug);
                xs = x_vec(is);
                ys = y_vec(js);
                xyz = [xs, ys];
            case 3
                [ind, is, js, ks] = tolStar(bli_tolerance, kgrid, point, debug);
                xs = x_vec(is);
                ys = y_vec(js);
                zs = z_vec(ks);
                xyz = [xs, ys, zs];
        end

        if mask_only
            
            % add current points to the mask
            mask(ind) = true;
            
        else
            
            % evaluate a BLI centered on point at grid nodes XYZ
            if scalar_dxyz
                if single_precision
                    mask_t = sinc( single(pi_on_dxyz .* bsxfun(@minus, xyz, point.')) );
                else
                    mask_t = sinc( pi_on_dxyz .* bsxfun(@minus, xyz, point.') );
                end
            else
                if single_precision
                    mask_t = sinc( single(bsxfun(@times, pi_on_dxyz, bsxfun(@minus, xyz, point.'))) );
                else
                    mask_t = sinc( bsxfun(@times, pi_on_dxyz, bsxfun(@minus, xyz, point.')) );
                end
            end
            mask_t = prod(mask_t, 2);

            % apply scaling for non-uniform grid
            if kgrid.nonuniform
                mask_t = mask_t * BLIscale;
            end

            % add this contribution to the overall source mask
            mask(ind) = mask(ind) + scale(point_ind) .* mask_t;
            
        end
        
    end
    
    % update the waitbar
    if display_wait_bar && (rem(point_ind, wait_bar_update_freq) == 0)
        waitbar(point_ind/num_points, h)
    end
    
end

% close the wait bar
if display_wait_bar
    close(h);
end