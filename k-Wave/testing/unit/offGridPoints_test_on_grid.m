function test_pass = offGridPoints_test_on_grid(plot_comparisons, ~)
%OFFGRIDPOINTS_TEST_ON_GRID Check truncated sinc approximation when points lie on the grid.
%
% DESCRIPTION:
%     offGridPoints_test_on_grid creates an off-grid point source using 
%     a sinc approximation and checks that when the point lies on the grid,
%     the off-grid mask does not extend.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 13th July 2021
%     last update      - 31st July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2021- Bradley Treeby

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

% set pass variable
test_pass = true;

% =========================================================================
% 1D - EVEN
% =========================================================================

% setup grid
Nx = 100;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% on grid point
scale = 1;
points_ind = Nx/2 - 1;
points = kgrid.x_vec(points_ind);
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass = false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(points_ind) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    plot(mask, '.');
    hold on;
    mask_tar = zeros(size(mask));
    mask_tar(points_ind) = 1;
    plot(mask_tar, 'ro');
    title('Mask');
    subplot(2, 1, 2);
    plot(pts);
    title('Values');
    sgtitle('On Grid Point');
end

% -------

% off grid point
scale = 1;
points = kgrid.x_vec(Nx/2 - 1) + kgrid.dx/4;
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    plot(mask, '.');
    title('Mask');
    subplot(2, 1, 2);
    plot(pts);
    title('Values');
    sgtitle('Off Grid Point');
end

% =========================================================================
% 1D - ODD
% =========================================================================

% setup grid
Nx = 101;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx);

% on grid point
scale = 1;
points_ind = ceil(Nx/2);
points = kgrid.x_vec(points_ind);
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(points_ind) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    plot(mask, '.');
    hold on;
    mask_tar = zeros(size(mask));
    mask_tar(points_ind) = 1;
    plot(mask_tar, 'ro');    
    title('Mask');
    subplot(2, 1, 2);
    plot(pts);
    title('Values');
    sgtitle('On Grid Point');
end

% -------

% off grid point
scale = 1;
points = kgrid.x_vec(ceil(Nx/2)) + kgrid.dx/4;
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(2, 1, 1);
    plot(mask, '.');
    title('Mask');
    subplot(2, 1, 2);
    plot(pts);
    title('Values');
    sgtitle('Off Grid Point');
end

% =========================================================================
% 2D - EVEN
% =========================================================================

% setup grid
Nx = 100;
Ny = 80;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% on grid point
scale = 1;
points_ind = [Nx/2 - 1, Ny/2 - 3];
points = [kgrid.x_vec(points_ind(1)), kgrid.y_vec(points_ind(2))].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(points_ind(1), points_ind(2)) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values');
    sgtitle('On Grid in X, On Grid in Y');
end

% -------

% off grid in x, on grid in y
scale = 1;
points = [kgrid.x_vec(Nx/2 - 1) + kgrid.dx/4, kgrid.y_vec(Ny/2 - 3)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('Off Grid in X, On Grid in Y');
end

% -------

% on grid in x, off grid in y
scale = 1;
points = [kgrid.x_vec(Nx/2 - 1), kgrid.y_vec(Ny/2 - 3) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('On Grid in X, Off Grid in Y');
end

% -------

% off grid in x, off grid in y
scale = 1;
points = [kgrid.x_vec(Nx/2 - 1) + kgrid.dx/4, kgrid.y_vec(Ny/2 - 3) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('Off Grid in X, Off Grid in Y');
end

% =========================================================================
% 2D - ODD
% =========================================================================

% setup grid
Nx = 101;
Ny = 81;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% on grid point
scale = 1;
points_ind = [ceil(Nx/2) - 1, ceil(Ny/2) - 3];
points = [kgrid.x_vec(points_ind(1)), kgrid.y_vec(points_ind(2))].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(points_ind(1), points_ind(2)) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values');
    sgtitle('On Grid in X, On Grid in Y');
end

% -------

% off grid in x, on grid in y
scale = 1;
points = [kgrid.x_vec(ceil(Nx/2) - 1) + kgrid.dx/4, kgrid.y_vec(ceil(Ny/2) - 3)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('Off Grid in X, On Grid in Y');
end

% -------

% on grid in x, off grid in y
scale = 1;
points = [kgrid.x_vec(ceil(Nx/2) - 1), kgrid.y_vec(ceil(Ny/2) - 3) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('On Grid in X, Off Grid in Y');
end

% -------

% off grid in x, off grid in y
scale = 1;
points = [kgrid.x_vec(ceil(Nx/2) - 1) + kgrid.dx/4, kgrid.y_vec(ceil(Ny/2) - 3) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    subplot(1, 2, 1);
    imagesc(mask);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Mask');
    
    subplot(1, 2, 2);
    imagesc(pts);
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');    
    title('Values');
    sgtitle('Off Grid in X, Off Grid in Y');
end

% =========================================================================
% 3D - EVEN
% =========================================================================

% setup grid
Nx = 100;
Ny = 80;
Nz = 60;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% on grid point
scale = 1;
x_ind = ceil(Nx/2) - 1;
y_ind = ceil(Ny/2) - 3;
z_ind = ceil(Nz/2) - 4;

points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind), kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(x_ind, y_ind, z_ind) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, On Grid in Y, On Grid in Z');
end

% -------

% off grid in x, on grid in y, on grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind), kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, On Grid in Y, On Grid in Z');
end

% -------

% on grid in x, off grid in y, on grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, Off Grid in Y, On Grid in Z');
end

% -------

% on grid in x, on grid in y, off grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind), kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, On Grid in Y, Off Grid in Z');
end

% -------

% off grid in x, off grid in y, on grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, Off Grid in Y, On Grid in Z');
end

% -------

% off grid in x, on grid in y, off grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind), kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, On Grid in Y, Off Grid in Z');
end

% -------

% on grid in x, off grid in y, off grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, Off Grid in Y, Off Grid in Z');
end

% -------

% off grid in x, off grid in y, off grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, Off Grid in Y, Off Grid in Z');
end

% =========================================================================
% 3D - ODD
% =========================================================================

% setup grid
Nx = 101;
Ny = 81;
Nz = 61;
dx = 1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% on grid point
scale = 1;
x_ind = ceil(Nx/2) - 1;
y_ind = ceil(Ny/2) - 3;
z_ind = ceil(Nz/2) - 4;

points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind), kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there is only one source point in the mask
if sum(mask(:)) > 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% check that the pts has a value of 1
if round(max(pts(:)), 5) ~= 1
    test_pass = false;
end

% check that the point is at the required position
if mask(x_ind, y_ind, z_ind) ~= 1
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, On Grid in Y, On Grid in Z');
end

% -------

% off grid in x, on grid in y, on grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind), kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, On Grid in Y, On Grid in Z');
end

% -------

% on grid in x, off grid in y, on grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, Off Grid in Y, On Grid in Z');
end

% -------

% on grid in x, on grid in y, off grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind), kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, On Grid in Y, Off Grid in Z');
end

% -------

% off grid in x, off grid in y, on grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind)].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, Off Grid in Y, On Grid in Z');
end

% -------

% off grid in x, on grid in y, off grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind), kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, On Grid in Y, Off Grid in Z');
end

% -------

% on grid in x, off grid in y, off grid in z
points = [kgrid.x_vec(x_ind), kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('On Grid in X, Off Grid in Y, Off Grid in Z');
end

% -------

% off grid in x, off grid in y, off grid in z
points = [kgrid.x_vec(x_ind) + kgrid.dx/4, kgrid.y_vec(y_ind) + kgrid.dx/4, kgrid.z_vec(z_ind) + kgrid.dx/4].';
mask = offGridPoints(kgrid, points, scale, 'MaskOnly', true, 'Debug', plot_comparisons);
pts  = offGridPoints(kgrid, points, scale, 'Debug', plot_comparisons);

% check there more than source point in the mask
if sum(mask(:)) == 1
    test_pass = false;
end

% check that mask and pts have the same number of entries
if sum(mask ~= 0) ~= sum(pts ~= 0)
    test_pass_false;
end

% plot
if plot_comparisons
    figure;
    
    subplot(2, 3, 1);
    imagesc(mask(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Mask - XY');
    
    subplot(2, 3, 2);
    imagesc(squeeze(mask(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Mask - XZ');    
    
    subplot(2, 3, 3);
    imagesc(squeeze(mask(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Mask - YZ');
    
    subplot(2, 3, 4);
    imagesc(pts(:, :, z_ind));
    axis image;
    xlabel('Y [grid points]');
    ylabel('X [grid points]');
    title('Values - XY');
    
    subplot(2, 3, 5);
    imagesc(squeeze(pts(:, y_ind, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('X [grid points]');
    title('Values - XZ');    
    
    subplot(2, 3, 6);
    imagesc(squeeze(pts(x_ind, :, :)));
    axis image;
    xlabel('Z [grid points]');
    ylabel('Y [grid points]');
    title('Values - YZ');
    
    sgtitle('Off Grid in X, Off Grid in Y, Off Grid in Z');
end