function test_pass = angularSpectrumCW_projection(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to propagate a piston source forwards and backwards.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th January 2019
%     last update - 19th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set comparison threshold
COMPARISON_THRESH = 0.5;

% set pass variable
test_pass = true;

% simulation settings
Nx = 128;
Ny = 128;
Nz = 32;
dx = 1e-4;
f0 = 2e6;
c0 = 1500;

% generate input plane as a Gaussian blob at the origin
FWHM = 15*dx;
variance = (FWHM / 2)^2 / 2 * log(2);
mean = 0;
kgrid = kWaveGrid(Nx, dx, Ny, dx);
r = sqrt(kgrid.x.^2 + kgrid.y.^2);
input_plane = exp(-(r - mean).^2 ./ (2 * variance));

% run forward simulation
output_plane = angularSpectrumCW(input_plane, dx, (Nz - 1) * dx, f0, c0);

% run reverse simulation
input_plane_remapped_1 = angularSpectrumCW(output_plane, dx, (Nz - 1) * dx, f0, c0, 'Reverse', true);
input_plane_remapped_2 = angularSpectrumCW(output_plane, dx, -(Nz - 1) * dx, f0, c0);

% compute relative L_inf error [%]
err_1 = 100 * max(abs(input_plane(:) - abs(input_plane_remapped_1(:)))) / max(abs(input_plane(:)));
err_2 = 100 * max(abs(input_plane(:) - abs(input_plane_remapped_2(:)))) / max(abs(input_plane(:)));

% check for test pass
if (err_1 > COMPARISON_THRESH) || (err_2 > COMPARISON_THRESH)
    test_pass = false;
end

% plot comparison
if plot_comparisons
    
    figure;
    subplot(4, 2, 1);
    imagesc(abs(input_plane));
    axis image;
    title('Input Plane - Amp');
    colorbar;
    
    subplot(4, 2, 2);
    imagesc(angle(input_plane));
    axis image;
    title('Input Plane - Phase');
    colorbar;
    
    subplot(4, 2, 3);
    imagesc(abs(output_plane));
    axis image;
    title('Output Plane - Amp');
    colorbar;    
    
    subplot(4, 2, 4);
    imagesc(angle(output_plane));
    axis image;
    title('Output Plane - Phase');
    colorbar;    
    
    subplot(4, 2, 5);
    imagesc(abs(input_plane_remapped_1));
    axis image;
    title('Re-projected Input Plane - Amp');
    colorbar;
    
    subplot(4, 2, 6);
    imagesc(angle(input_plane_remapped_1));
    axis image;
    title('Re-projected Input Plane - Phase');
    colorbar;
    
    subplot(4, 2, 7);
    imagesc(100 * abs(abs(input_plane) - abs(input_plane_remapped_1)) ./ max(abs(input_plane(:))) );
    axis image;
    title('Error - Amp [%]');
    colorbar;
    
    subplot(4, 2, 8);
    imagesc(abs(angle(input_plane) - angle(input_plane_remapped_1)));
    axis image;
    title('Error - Phase [deg]');
    colorbar;
    
    scaleFig(1, 1.5);
    
end