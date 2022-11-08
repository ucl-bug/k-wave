function affine = getAffineMatrix(translation, rotation)
%GETAFFINEMATRIX Return matrix for affine transform in 3D.
%
% DESCRIPTION:
%     getAffineMatrix returns an affine matrix defined by the specified
%     translation and rotation. 
%
% USAGE:
%     affine = getAffineMatrix(translation, rotation)
%
% INPUTS:
%     translation         - translation given as [dx, dy] in 2D and 
%                           [dx, dy, dz] in 3D
%     rotation            - rotation angle/s in degrees given as 
%                           [th] in 2D (counter-clockwise) and 
%                           [x_th, y_th, z_th] in 3D (rotate about x
%                           then y' then z'') [degrees] 
%
% OUTPUTS:
%     affine              - affine transformation matrix
%
% ABOUT:
%     author              - Bradley Treeby
%     date                - 13th September 2018
%     last update         - 13th September 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Bradley Treeby
%
% See also computeLinearTransform

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

% check dimensions
if (numel(translation) == 2) && (numel(rotation) == 1)
    
    % assign the inputs
    dx    = translation(1);
    dy    = translation(2);
    th    = rotation;
    
    % build affine matrix (counter-clockwise)
    affine = [ cosd(th) -sind(th)  dx
               sind(th)  cosd(th)  dy
                    0         0    1];
    
elseif (numel(translation) == 3) && (numel(rotation) == 3)
    
    % assign the inputs
    dx   = translation(1);
    dy   = translation(2);
    dz   = translation(3);
    x_th = rotation(1);
    y_th = rotation(2);
    z_th = rotation(3);

    % build the rotation matrices
    x_th_matrix = [ 1       0           0
                     0 cosd(x_th) -sind(x_th)
                     0 sind(x_th)  cosd(x_th)];

    y_th_matrix = [ cosd(y_th) 0 sind(y_th)
                         0     1      0
                   -sind(y_th) 0 cosd(y_th)];

    z_th_matrix = [ cosd(z_th) -sind(z_th) 0
                    sind(z_th)  cosd(z_th) 0
                         0           0     1];

    % build affine matrix
    affine = zeros(4, 4);
    affine(1:3, 1:3) = z_th_matrix * y_th_matrix * x_th_matrix;
    affine(:, 4) = [dx, dy, dz, 1];
    
else
    error('Incorrect size for translation and rotation inputs.');
end