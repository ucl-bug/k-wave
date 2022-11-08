function [bowls, bowls_labelled] = makeMultiBowl(grid_size, bowl_pos, radius, diameter, focus_pos, varargin)
%MAKEMULTIBOWL Create a binary map of multiple bowls within a 3D grid.
%
% DESCRIPTION:
%     makeMultiBowl creates a binary map of multiple bowls within a
%     three-dimensional grid using makeBowl. The position of the bowls is
%     denoted by 1's in the the matrix with 0's elsewhere. A labelled
%     matrix can also be returned, where the position of the first bowl
%     is denoted by 1's, the position of the second bowl by 2's, and so
%     on.
%
% USAGE:
%     [bowls, bowl_labelled] = makeMultiBowl(grid_size, bowl_pos, radius, diameter, focus_pos)
%     [bowls, bowl_labelled] = makeMultiBowl(grid_size, bowl_pos, radius, diameter, focus_pos, ...)
%
% INPUTS:
%     grid_size       - size of the 3D grid given as a three element
%                       vector [Nx, Ny, Nz] [grid points] 
%     bowl_pos        - centre of the rear surface of each bowl given as
%                       matrix of dimensions N x 3, with each row 
%                       specifying the centre for each bowl as a three
%                       element vector [bx, by, bz] [grid points]
%     radius          - radius of curvature of each bowl given as either
%                       a single number (if the bowls have the same
%                       radius), or an N-element vector containing the
%                       radius for each bowl [grid points]
%     diameter        - aperture diameter of each bowl given as either a
%                       single number (if the bowls have the same
%                       diameter), or an N-element vector containing the
%                       diameter for each bowl [grid points] 
%     focus_pos       - any point on the beam axis of the bowl given as
%                       either a three element vector (if the bowls have
%                       the same focus_pos), or as a matrix of dimensions 
%                       N x 3, with each row specifying the focus_pos for
%                       each bowl given as a three element vector
%                       [fx, fy, fz] [grid points]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Binary'        - Boolean controlling whether the bowl map is
%                       returned as a double precision matrix (false) or
%                       a logical matrix (true) (default = false).
%     'Plot'          - Boolean controlling whether the bowl array is
%                       plotted using voxelPlot (default = false).
%     'RemoveOverlap' - Boolean controlling whether overlapped grid points
%                       within each bowl (not between bowls) are removed
%                       (default = false). 
%
% OUTPUTS:
%     bowls           - 3D binary map of bowls
%     bowls_labelled  - 3D labelled matrix of bowls
%
% ABOUT:
%     author          - Yan To Ling and Bradley Treeby
%     date            - 17th November 2014
%     last update     - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Yan To Ling and Bradley Treeby
%
% See also makeBowl, makeMultiArc

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

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% optional input defaults
NUM_REQUIRED_INPUTS     = 5;
PLOT_BOWL_DEF           = false;
BINARY_DEF              = false;
REMOVE_OVERLAP_DEF      = false;

% =========================================================================
% INPUT CHECKING
% =========================================================================

% assign optional inputs with defaults
plot_bowl               = PLOT_BOWL_DEF;
binary                  = BINARY_DEF;
remove_overlap          = REMOVE_OVERLAP_DEF;

% replace with user defined values if provided
if nargin < NUM_REQUIRED_INPUTS
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Binary'
                
                % assign input
                binary = varargin{input_index + 1};
                
                % check value
                if ~islogical(binary)
                    error('Optional input ''Binary'' must be Boolean.');
                end
                
            case 'Plot'
                
                % assign input
                plot_bowl = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_bowl)
                    error('Optional input ''Plot'' must be Boolean.');
                end
                
            case 'RemoveOverlap'
                
                % assign input
                remove_overlap = varargin{input_index + 1};
                
                % check value
                if ~islogical(remove_overlap)
                    error('Optional input ''RemoveOverlap'' must be Boolean.');
                end
                
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check inputs
if size(bowl_pos, 2) ~= 3
    error('bowl_pos should contain 3 columns, with [bx, by, bz] in each row.');
end
if length(radius) ~= 1 && length(radius) ~= size(bowl_pos,1)
    error('The number of rows in bowl_pos and radius does not match.');
end
if length(diameter) ~= 1 && length(diameter) ~= size(bowl_pos,1)
    error('The number of rows in bowl_pos and diameter does not match.');
end

% force integer grid size values
grid_size = round(grid_size);
bowl_pos  = round(bowl_pos);
focus_pos = round(focus_pos);
diameter  = round(diameter);
radius    = round(radius);

% =========================================================================
% CREATE BOWLS
% =========================================================================

% preallocate output matrices
if binary
    bowls = false(grid_size);
else
    bowls = zeros(grid_size);
end
if nargout == 2
    bowls_labelled = zeros(grid_size);
end

% loop for calling makeBowl
for bowl_index = 1:size(bowl_pos, 1)
    
    % update command line status
    if bowl_index == 1
        tic;
    else
        toc; tic;
    end
    fprintf(['Creating bowl ' num2str(bowl_index) ' of ' num2str(size(bowl_pos,1)) '... ']);
    
    % get parameters for current bowl
    if size(bowl_pos, 1) > 1
        bowl_pos_k = bowl_pos(bowl_index, :);
    else
        bowl_pos_k = bowl_pos;
    end
    if length(radius) > 1
        radius_k = radius(bowl_index);
    else
        radius_k = radius;
    end
    if length(diameter) > 1
        diameter_k = diameter(bowl_index);
    else
        diameter_k = diameter;
    end      
    if size(focus_pos, 1) > 1
        focus_pos_k = focus_pos(bowl_index, :);
    else
        focus_pos_k = focus_pos;
    end    
    
    % create new bowl
    new_bowl = makeBowl(grid_size, bowl_pos_k, radius_k,...
        diameter_k, focus_pos_k, 'RemoveOverlap', remove_overlap, ...
        'Binary', binary);
        
    % add bowl to bowl matrix
    bowls = bowls + new_bowl;
    
    % add new bowl to labelling matrix
    if nargout == 2
        bowls_labelled(new_bowl == 1) = bowl_index;
    end
    
end

toc;

% check if any of the bowls are overlapping
if maxND(bowls) > 1
    
    % display warning
    disp(['WARNING: ' num2str(maxND(bowls) - 1) ' bowls are overlapping']);
    
    % force the output to be binary
    bowls(bowls ~= 0) = 1;
    
end

% create the figure
if plot_bowl
    voxelPlot(bowls);
end