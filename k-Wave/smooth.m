function mat_sm = smooth(varargin)
%SMOOTH Smooth a matrix.
%
% DESCRIPTION:
%     smooth filters an input matrix using an n-dimensional frequency
%     domain window created using getWin. If no window type is specified, a
%     Blackman window is used. 
%
% USAGE:
%     mat_sm = smooth(mat)
%     mat_sm = smooth(mat, restore_max)
%     mat_sm = smooth(mat, [], window_type)
%     mat_sm = smooth(mat, restore_max, window_type)
%
% INPUTS:
%     mat         - spatial distribution to smooth
%
% OPTIONAL INPUTS:
%     restore_max - Boolean controlling whether the maximum value is
%                   restored after smoothing (default = false).
%     window_type - shape of the smoothing window; any valid inputs to
%                   getWin are supported (default = 'Blackman').
%
% OUTPUTS:
%     mat_sm      - smoothed spatial distribution
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 29th April 2009
%     last update - 29th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby
%
% See also fft, ifft, fft2, ifft2, fftn, ifftn, kWaveGrid, getWin, numDim

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

% define literals
DEF_WINDOW_TYPE = 'Blackman';
DEF_USE_ROTATION = true;
DEF_RESTORE_MAX = false;

% check if first option is a kWaveGrid object
if isa(varargin{1}, 'kWaveGrid')
    
    % give a warning
    disp('Warning: inputs to smooth have changed (kgrid input is no longer required)');
    
    % call smooth recursively without kgrid input
    mat_sm = smooth(varargin{2:end});
    
    % exit
    return
    
end

% assign input
mat = varargin{1};

% check value
validateattributes(mat, {'numeric'}, {'finite'}, 'smooth', 'mat', 1);
   
% check for restore max input
if nargin < 2 || isempty(varargin{2})
    
    % assign default
    restore_max = DEF_RESTORE_MAX;
    
else
    
    % assign input
    restore_max = varargin{2};
    
    % check value
    validateattributes(restore_max, {'logical'}, {'scalar'}, 'smooth', 'restore_max', 2);
    
end

% check for window type input
if nargin < 3 || isempty(varargin{3})
    
    % assign default
    window_type = DEF_WINDOW_TYPE;
    
else
    
    % assign input
    window_type = varargin{3};
    
    % check value
    validateattributes(window_type, {'char'}, {}, 'smooth', 'window_type', 3);
    
end

% get the grid size
grid_size = size(mat);

% remove singleton dimensions
if numDim(mat) ~= length(grid_size)
    grid_size(grid_size == 1) = [];
end

% use a symmetric filter for odd grid sizes, and a non-symmetric filter for
% even grid sizes to ensure the DC component of the window has a value of
% unity
window_symmetry = logical(rem(grid_size, 2));

% get the window, taking the absolute value to discard machine precision
% negative values 
win = abs(getWin(grid_size, window_type, 'Rotation', DEF_USE_ROTATION, 'Symmetric', window_symmetry));

% rotate window if input mat is (1, N)
if isrow(mat)
    win = win.';
end

% apply the filter
mat_sm = real(ifftn(fftn(mat) .* ifftshift(win)));

% restore magnitude if required
if restore_max
    mat_sm = ( max(abs(mat(:))) ./ max(abs(mat_sm(:))) ) .* mat_sm;
end