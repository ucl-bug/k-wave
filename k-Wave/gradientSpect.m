function varargout = gradientSpect(f, dn, dim, deriv_order)
%GRADIENTSPECT Calculate the gradient using a Fourier spectral method.
%
% DESCRIPTION:
%     gradientSpect calculates the gradient of an n-dimensional input
%     matrix using the Fourier collocation spectral method. The gradient
%     for singleton dimensions is returned as 0.
%
%     Examples:
%     1. 1D Sinusoid:
%
%         % compute gradient of a 2 period sinusoid
%         x = pi/20:pi/20:4*pi;
%         y = sin(x);
%         dydx = gradientSpect(y, pi/20);
%
%         % plot gradient and error compared to analytical solution
%         subplot(2, 1, 1);
%         plot(x, cos(x), 'k-', x, dydx, 'bx');
%         axis tight;
%         title('dy/dx');
%         subplot(2, 1, 2);
%         plot(x, cos(x) - dydx, 'k-');
%         axis tight;
%         title('Relative Error');
%
%     2. Modified 2D Mathworks gradient example (x and y reversed):
%
%         % generate input and compute gradient
%         [x, y] = meshgrid(-2:.2:2, -2:.2:2);
%         z = x .* exp(-x.^2 - y.^2);
%         [px, py] = gradientSpect(z, [.2, .2]);
%
%         % plot
%         figure;
%         contour(z);
%         hold on;
%         quiver(py, px);
%         hold off;
%
% USAGE:
%     fx = gradientSpect(f, dx)
%     fx = gradientSpect(f, dx, [], deriv_order)
%     fn = gradientSpect(f, dn, dim)
%     fn = gradientSpect(f, dn, dim, deriv_order)
%     [fx, fy] = gradientSpect(f, dn)
%     [fx, fy] = gradientSpect(f, dn, [], deriv_order)
%     [fx, fy, fz, ...] = gradientSpect(f, dn)
%     [fx, fy, fz, ...] = gradientSpect(f, dn, [], deriv_order)
%
% INPUTS:
%     f           - matrix or vector to find the gradient of
%     dn          - array of values for the grid point spacing in each
%                   dimension. If a value for dim is given, dn is the
%                   spacing in dimension dim.
%       
% OPTIONAL INPUTS:
%     dim         - optional input to specify a single dimension over which
%                   to compute the gradient for n-dimension input functions
%     deriv_order - order of the derivative to compute, e.g., use 1 to
%                   compute df/dx, 2 to compute df^2/dx^2, etc 
%                   (default = 1)
%
% OUTPUTS:
%     fx, fy, ... - gradient in the each dimension, where x corresponds
%                   to dim = 1, y corresponds to dim = 2 etc 
%       
% ABOUT:
%     author      - Bradley Treeby
%     date        - 25th January 2012
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby
%
% See also gradient, gradientFD

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

%#ok<*AGROW>

% check for user input for deriv_order
if nargin < 4  || isempty(deriv_order)
    deriv_order = 1;
elseif ~((deriv_order > 0) && (deriv_order == round(deriv_order)))
    error('input for deriv_order must be an integer > 0.');
end

% check for user input for dim
if (nargin < 3) || isempty(dim)
    dim = 0;
end

% get size of the input function
sz = size(f);

% check if input is 1D or user defined input dimension is given
if dim || (nargout == 1 && (max(sz) == prod(sz)))
    
    % check if a single dn value is given, if not, extract the required
    % value
    if numel(dn) ~= 1
        dn = dn(dim);
    end
    
    % get the grid size along the specified dimension, or the longest
    % dimension if 1D
    if nargout == 1 && (max(sz) == prod(sz))
        [Nx, dim] = max(sz);
    else
        Nx = sz(dim);
    end
    
    % get the wavenumbers
    kx = getWavenumbers(Nx, dn, dim);

    % calculate derivative and assign output       
    varargout{1} = ifft( bsxfun(@times, (1i * kx).^deriv_order, fft(f, [], dim)) , [], dim, 'symmetric');
    
else
    
    % check for the correct number of dn values
    if length(dn) ~= length(sz)
        error([num2str(length(sz)) ' values for dn must be specified for a ' num2str(length(sz)) '-dimensional input matrix.']);
    end
    
    % check for the correct number of outputs
    if nargout ~= length(sz)
        error(['Incorrect number of output arguments for ' num2str(length(sz)) '-dimensional input matrix.']);
    end

    % calculate the gradient for each non-singleton dimension
    for dim = 1:ndims(f)
        if sz(dim) > 1

            % get the wavenumbers
            kx = getWavenumbers(sz(dim), dn(dim), dim);

            % calculate derivative and assign output       
            varargout{dim} = ifft( bsxfun(@times, (1i * kx).^deriv_order, fft(f, [], dim)) , [], dim, 'symmetric'); 

        else
            
            % assign the derivate for singleton dimensions to be 0
            varargout{dim} = 0;
            
        end
    end
end

function kx = getWavenumbers(Nx, dx, dim)
% Subfuction to compute the wavenumber vector.

% calculate the wavenumbers
if rem(Nx, 2) == 0
    
    % grid dimension has an even number of points
    nx = ((-Nx/2:Nx/2-1)/Nx).';
    
else
    
    % grid dimension has an odd number of points
    nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx).';
    
end
kx = ifftshift((2*pi/dx).*nx); 

% permute to be in the correction direction for use with bsxfun
if dim == 1
    kx = reshape(kx, Nx, 1);
else
    kx = reshape(kx, [ones(1, dim - 1), Nx]);
end
