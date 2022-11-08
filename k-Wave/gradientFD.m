function varargout = gradientFD(f, dn, dim, deriv_order, accuracy_order)
%GRADIENTFD Calculate the gradient using a finite-difference method.
%
% DESCRIPTION:
%     gradientFD calculates the gradient of an n-dimensional input matrix
%     using the finite-difference method. For one-dimensional inputs, the
%     gradient is always computed along the non-singleton dimension. For
%     higher dimensional inputs, the gradient for singleton dimensions is
%     returned as 0. For elements in the centre of the grid, the gradient
%     is computed using centered finite-differences. For elements on the
%     edge of the grid, the gradient is computed using forward or backward
%     finite-differences. The order of accuracy of the finite-difference
%     approximation is controlled by accuracy_order (default = 2). The
%     calculations are done using sparse multiplication, so the input
%     matrix is always cast to double precision.  
%
% USAGE:
%     fx = gradientFD(f, dx)
%     fx = gradientFD(f, dx, [], deriv_order)
%     fx = gradientFD(f, dx, [], deriv_order, accuracy_order)
%     fn = gradientFD(f, dn, dim)
%     fn = gradientFD(f, dn, dim, deriv_order, accuracy_order)
%     [fx, fy] = gradientFD(f, dn)
%     [fx, fy] = gradientFD(f, dn, [], deriv_order, accuracy_order)
%     [fx, fy, fz, ...] = gradientFD(f, dn)
%     [fx, fy, fz, ...] = gradientFD(f, dn, [], deriv_order, accuracy_order)
%
% INPUTS:
%     f               - matrix or vector to find the gradient of
%     dn              - array of values for the grid point spacing in each
%                       dimension. If a value for dim is given, dn is the
%                       spacing in dimension dim.
%       
% OPTIONAL INPUTS:
%     dim             - optional input to specify a single dimension over
%                       which to compute the gradient for n-dimension
%                       input functions
%     deriv_order     - order of the derivative to compute, e.g., use 1 to
%                       compute df/dx, 2 to compute df^2/dx^2, etc.
%                       (default = 1)
%     accuracy_order  - order of accuracy for the finite difference
%                       coefficients. Because centered differences are
%                       used, this must be set to an integer multiple of
%                       2 (default = 2)
%       
% OUTPUTS:
%     fx, fy, ...     - gradient in the each dimension, where x corresponds
%                       to dim = 1, y corresponds to dim = 2 etc 
%       
% ABOUT:
%       author        - Bradley Treeby
%       date          - 24th August 2012
%       last update   - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby
%
% See also getFDMatrix, gradient, gradientSpect

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

%#ok<*ISMAT>
%#ok<*AGROW>

% % display warning if input is not double precision
% if ~isa(f, 'double')
%     disp('gradientFD: converting input to double precision');
% end

% force input to be in double precision (sparse operations are only
% supported in double or logical formats)
f = double(f);

% get size of the input function
sz = size(f);

% check dimension of input
f_dims = numDim(f);

% check for user input for accuracy_order
if nargin < 5 || isempty(accuracy_order)
    accuracy_order = 2;
elseif rem(accuracy_order, 2)
    error('Input for accuracy_order must be an integer multiple of 2.');
end

% check for user input for deriv_order
if nargin < 4 || isempty(deriv_order)
    deriv_order = 1;
elseif ~((deriv_order > 0) && (deriv_order == round(deriv_order)))
    error('Input for deriv_order must be an integer > 0.');
end

% check for user input for dim
if (nargin < 3) || isempty(dim)
    
    % check if input is 1D
    if f_dims == 1
        
        % only compute gradient over required dimension
        if sz(1) == 1
            dim_array = 2;
        else
            dim_array = 1;
        end
        
    else
        
        % otherwise compute gradient over all dimensions
        dim = 0;
        dim_array = 1:f_dims;
    
        % check for the correct number of dn values
        if length(dn) ~= length(sz)
            error([num2str(length(sz)) ' values for dn must be specified for a ' num2str(length(sz)) '-dimensional input matrix.']);
        end
        
    end
else
    dim_array = dim;
end

% only allow 1, 2, or 3D inputs (only these cases are implemented)
if f_dims > 3
    error('Input for f must have 1, 2, or 3 dimensions.');
end

% if input is 1D, only allow it to be a row or column vector
if (f_dims == 1) && (length(sz) > 2)
    error('1D inputs for f must be row or column vectors.');
end

% if input is 2D, only allow it to be a regular matrix
if (f_dims == 2) && (length(sz) > 2)
    error('2D inputs for f must be of size m x n.');
end

% if input is 2D or 3D, check dim input isn't bigger than the matrix size
if (f_dims > 1) && (dim > f_dims)
    error('Input for dim cannot be greater than the number of dimensions of f.');
end

% set output argument index
argout_index = 1;

% loop through required dimensions
for dim_index = dim_array

    % set dimension
    dim = dim_index;
    
    % check if a single dn value is given, if not, extract the required
    % value
    if numel(dn) ~= 1
        dn_val = dn(dim);
    else
        dn_val = dn;
    end
    
    % get the grid size along the specified dimension, or the longest
    % dimension if 1D
    if nargout == 1 && (max(sz) == prod(sz))
        [Nx, dim] = max(sz);
    else
        Nx = sz(dim);
    end
    
    % get the finite-difference matrix
    FDM = getFDMatrix(Nx, dn_val, deriv_order, accuracy_order);
    
    % compute derivatives of 1D or 2D matrix
    if ndims(f) < 3
        switch dim

            % derivative along dimension 1 (columns)
            case 1
                varargout{argout_index} = FDM * f; 

            % derivative along dimension 2 (rows)
            case 2
                varargout{argout_index} = f * FDM.';

        end

    % compute derivatives of 3D matrix    
    else
        switch dim

            % derivative along dimension 1
            case 1

                % preallocate output matrix
                varargout{argout_index} = zeros(size(f));

                % loop through z-plane
                for dim3_index = 1:sz(3)
                    varargout{argout_index}(:, :, dim3_index) = FDM * f(:, :, dim3_index);
                end

            % derivative along dimension 2
            case 2

                % preallocate output matrix
                varargout{argout_index} = zeros(size(f));

                % rotate FDM
                FDM = FDM.';

                % loop through z-plane
                for dim3_index = 1:sz(3)
                    varargout{argout_index}(:, :, dim3_index) = f(:, :, dim3_index) * FDM;
                end

            % derivative along dimension 3   
            case 3

                % preallocate output matrix
                varargout{argout_index} = zeros(size(f));

                % rotate FDM
                FDM = FDM.';

                % loop through x-plane
                for dim1_index = 1:sz(1)
                    varargout{argout_index}(dim1_index, :, :) = squeeze(f(dim1_index, :, :)) * FDM;
                end

        end
        
    end
            
    % increment output argument index
    argout_index = argout_index + 1;
    
end