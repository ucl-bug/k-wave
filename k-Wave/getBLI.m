function varargout = getBLI(func, dx, upsampling_factor, plot_BLI)
%GETBLI Compute underlying Fourier band-limited interpolant (BLI).
%
% DESCRIPTION:
%     getBLI computes the underlying band-limited interpolant (BLI) of a
%     discretely sampled 1D function assuming a Fourier basis.
%
%     Example: Plot the BLI of a discrete delta function
%
%         delta = zeros(20, 1);
%         delta(end/2) = 1;
%         getBLI(delta, 0.001, 10, true);
%
% USAGE:
%     BLI = getBLI(func)
%     BLI = getBLI(func, dx)
%     BLI = getBLI(func, dx, sampling_factor)
%     BLI = getBLI(func, dx, sampling_factor, plot_BLI)
%     BLI = getBLI(func, [], [], plot_BLI)
%     [x_fine, BLI] = getBLI(func, dx, sampling_factor, plot_BLI)
%
% INPUTS:
%     func              - 1D input function
%
% OPTIONAL INPUTS:
%     dx                - spatial sampling [m] (default = 1)
%     upsampling_factor - upsampling factor used to sample the underlying
%                         BLI (default = 20)
%     plot_BLI          - Boolean controlling whether the discrete input
%                         function and its BLI are displayed
%                         (default = false)   
%
% OUTPUTS:
%     BLI               - band-limited interpolant
%     x_fine            - x-grid for BLI 
%
% ABOUT:
%     author            - Bradley Treeby
%     date              - 4th July 2014
%     last update       - 13th January 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2019 Bradley Treeby
%
% See also getDeltaBLI, spectrum

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

% set defaults
DX_DEF = 1;
UPSAMPLING_FACTOR_DEF = 20;
PLOT_BLI_DEF = false;

% check inputs
if nargin < 2 || isempty(dx)
    dx = DX_DEF;
end
if nargin < 3 || isempty(upsampling_factor)
    upsampling_factor = UPSAMPLING_FACTOR_DEF;
end
if nargin < 4 || isempty(plot_BLI)
    plot_BLI = PLOT_BLI_DEF;
end

% force input function to be in desired direction
func = reshape(func, 1, []);

% extract parameters
Nx = length(func);

% compute set of wavenumbers
dk = 2 * pi / (dx * Nx);
if rem(Nx, 2)
    
    % Nx is odd 
    k_min = -pi / dx + dk / 2;
    k_max =  pi / dx - dk / 2;
    
else
    
    % Nx is even
    k_min = -pi / dx;
    k_max =  pi / dx - dk;
    
end
k = k_min:dk:k_max;

% compute fine grid
x_fine = 0:(dx / upsampling_factor):((Nx - 1) * dx);

% compute Fourier coefficients
func_k = fftshift(fft(func)) / Nx;

% form 
BLI = real(sum( func_k * exp(1i * k.' * x_fine), 1));

% assign output
if nargout == 2
    varargout{1} = x_fine;
    varargout{2} = BLI;
else
    varargout{1} = BLI;
end

% plot BLI if desired
if nargin == 4 && plot_BLI
    
    % get suitable plot scale
    [x_sc, scale, prefix] = scaleSI((Nx - 1) * dx); %#ok<ASGLU>
    
    % create figure
    figure;
    plot(scale * (0:dx:((Nx - 1) * dx)), func, 'ro')
    hold on;
    plot(scale * x_fine, BLI, 'k-');
    
    % adjust plot scale
    axis tight;
    plot_range = max(BLI(:)) - min(BLI(:));
    set(gca, 'YLim', [min(BLI(:)) - 0.05 * plot_range, max(BLI(:)) + 0.05 * plot_range]);
    
    % annotate
    legend('Discrete Function', 'Band-Limited Interpolant', 'Location', 'Best');
    xlabel(['x [' prefix 'm]']);
    ylabel('Amplitude');
    
end