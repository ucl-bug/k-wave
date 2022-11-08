function voxelPlot(mat, varargin)
%VOXELPLOT 3D plot of voxels in a binary matrix.
%
% DESCRIPTION:
%     voxelPlot produces a 3D plot of a binary matrix, where filled voxels
%     are displayed at the positions of the 1's. The colormap,
%     transparency, and axis limits can be controlled through optional
%     inputs. The input matrix must be in single or double precision.
%
%     Examples:
%         voxelPlot(makeBall(30, 30, 30, 15, 15, 15, 12));
%         voxelPlot(makeBall(20, 20, 20, 10, 10, 10, 4), 'AxisTight', true, 'Color', [1 0 0], 'Transparency', 0.5);
%
% USAGE:
%     voxelPlot(mat)
%     voxelPlot(mat, ...)
%
% INPUTS:
%     mat         - input 3D matrix in single or double precision
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'AxisTight' - Boolean controlling whether axis limits are set to
%                   only display the filled voxels (default = false).
%     'Color'     - Three element array specifying rgb color (default =
%                   [1, 1, 0.4]).
%     'Transparency' 
%                 - value between 0 and 1 specifying transparency where 1
%                   gives no transparency (default = 0.8) 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 3rd September 2009
%     last update - 7th June 2017
%
%     voxelPlot calls the function image3Ddata by Kevin Moerman from
%     MATLAB central (available from http://www.mathworks.com/...
%     matlabcentral/fileexchange/24081-image3ddata). image3Ddata is
%     redistributed with k-Wave under the terms of the BSD license. 
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also patch

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

% check input matrix is 3D and single or double precision
if numDim(mat) ~= 3 || ~isfloat(mat)
    error('Input must be a 3D matrix in single or double precision.');
end

% set literals
num_req_input_variables = 1;
transparency = 0.8;
axis_tight = false;
color_map = [1, 1, 0.4];    % yellow

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif rem(nargin - num_req_input_variables, 2)
    error('Optional input parameters must be given as param, value pairs.');    
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'AxisTight'
                axis_tight = varargin{input_index + 1}; 
            case 'Color'
                color_map  = varargin{input_index + 1};                 
            case 'Transparency'
               transparency = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% scale to a max of 1
mat = mat ./ max(mat(:));

% create structure array containing coordinate and colour data for 3D image
[IMAGE_3D_DATA] = image3Ddata(mat);  

% threshold and select the voxels to display
voxel_num = (mat == 1);  
voxel_face_num = IMAGE_3D_DATA.voxel_patch_face_numbers(voxel_num, :);  
M_faces = IMAGE_3D_DATA.voxel_patch_faces(voxel_face_num, :);  
M_vertices = IMAGE_3D_DATA.corner_coordinates_columns_XYZ;  

% create a new figure with a white background
fig = figure;
set(fig, 'Color', [1, 1, 1]); 

% plot the voxels using patch
hp2 = patch('Faces', M_faces, 'Vertices', M_vertices, 'EdgeColor', ...
    'black', 'CData', IMAGE_3D_DATA.voxel_patch_CData(voxel_face_num,:), ...
    'FaceColor', 'flat');  

% set the tranparency
set(hp2, 'FaceAlpha', transparency);

% set the axes properties and colormap
view(45, 30); 
axis equal;
box on;
colormap(color_map); 
caxis([0, 1]); 
grid on;  

% add the axes labels
xlabel('y [voxels]');
ylabel('x [voxels]');
zlabel('z [voxels]');

% force the display to be the same size as mat
if ~axis_tight
    sz = size(mat);
    set(gca, 'XLim', [0.5, sz(2) + 0.5], ...
             'YLim', [0.5, sz(1) + 0.5], ...
             'ZLim', [0.5, sz(3) + 0.5]);
end