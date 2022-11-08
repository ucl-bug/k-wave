function b_mode = scanConversion(scan_lines, steering_angles, image_size, c0, dt, resolution)
%SCANCONVERSION Convert scan-lines in polar coordinates to a B-mode ultrasound image.
%
% DESCRIPTION:
%     scanConversion computes the remapping between a series of scan lines
%     in polar coordinates (i.e., taken at different steering angles) to a
%     B-mode ultrasound image in Cartesian coordinates suitable for
%     display. The remapping is performed using bilinear interpolation via
%     interp2.
%
% USAGE:
%     b_mode = scanConversion(scan_lines, steering_angles, image_size, c0, dt)
%     b_mode = scanConversion(scan_lines, steering_angles, image_size, c0, dt, resolution)
%
% INPUTS:
%     scan_lines      - matrix of scan lines indexed as (angle, time)
%     steering_angles - array of scanning angles [degrees]
%     image_size      - size of the desired output image [x, y] in metres
%     c0              - sound speed in the medium [m/s]
%     dt              - time step used in the acquisition of the scan
%                       lines [s]
%     resolution      - optional input to set the resolution of the
%                       output images in pixels (default = [256, 256])
%
% OUTPUTS:
%     b_mode          - the converted B-mode ultrasound image
%
% ABOUT:
%     author          - Bradley E. Treeby
%     date            - 23rd February 2011
%     last update     - 7th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby
%
% See also interp2

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% define literals
X_RESOLUTION_DEF = 256;     % [pixels]
Y_RESOLUTION_DEF = 256;     % [pixels]

% check for the resolution inputs
if nargin == 6
    x_resolution = resolution(1);
    y_resolution = resolution(end);
else
    x_resolution = X_RESOLUTION_DEF;
    y_resolution = Y_RESOLUTION_DEF;
end

% assign the inputs
x = image_size(1);
y = image_size(2);

% start the timer
tic;

% update command line status
disp('Computing ultrasound scan conversion...');

% extract a_line parameters
Nt = length(scan_lines(1, :));

% calculate radius variable based on the sound speed in the medium and the
% round trip distance
r = c0 * (1:Nt) * dt / 2;     % [m]

% create regular Cartesian grid to remap to
pos_vec_y_new = (0:1 / (y_resolution-1):1) .* y - y / 2;
pos_vec_x_new = (0:1 / (x_resolution-1):1) .* x;
[pos_mat_y_new, pos_mat_x_new] = ndgrid(pos_vec_y_new, pos_vec_x_new);

% convert new points to polar coordinates
[th_cart, r_cart] = cart2pol(pos_mat_x_new, pos_mat_y_new);

% interpolate using linear interpolation
b_mode = interp2(r, 2 * pi * steering_angles ./ 360, scan_lines, r_cart, th_cart, 'linear').';

% set any values outside of the interpolation range to be 0
% b_mode(isnan(b_mode)) = max(b_mode(:));
b_mode(isnan(b_mode)) = 0;

% update command line status
disp(['  completed in ' scaleTime(toc)]);             
