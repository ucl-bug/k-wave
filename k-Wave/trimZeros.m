function [data, ind] = trimZeros(data)
%TRIMZEROS Create a tight bounding box by removing zeros.
%
% DESCRIPTION:
%     trimZeros reduces the size of a matrix to create a tight bounding box
%     around the data by removing zeros. For example, in 2D, if the input
%     data is given by:
%
%         0 0 0 0 0 0
%         0 0 0 3 0 0
%         0 0 1 3 4 0
%         0 0 1 3 4 0
%         0 0 1 3 0 0
%         0 0 0 0 0 0
%
%     trimZeros will return:
%
%         0 3 0
%         1 3 4
%         1 3 4
%         1 3 0
%
% USAGE:
%     data = trimZeros(data)
%     [data, ind] = trimZeros(data)
%
% INPUTS:
%     data        - matrix to trim
%
% OUTPUTS:
%     data        - trimmed matrix
%     ind         - indices used to trim matrix
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 31st May 2017
%     last update - 9th December 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2022 Bradley Treeby
%
% See also expandMatrix, resize

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

% only allow 1D, 2D, and 3D
if numDim(data) > 3
    error('Input data must be 1D, 2D, or 3D.');
end

% set collapse directions for each dimension
collapse = [2, 3; 1, 3; 1, 2];

% preallocate output to store indices
ind = zeros(1, 2 * numDim(data));

% loop through dimensions
for dim_index = 1:numDim(data)

    % collapse to 1D vector
    if numDim(data) == 1
        summed_values = data;
    else
        summed_values = sum(sum(abs(data), collapse(dim_index, 1)), collapse(dim_index, 2));
    end
    
    % find the first and last non-empty values
    ind_first = find(summed_values > 0, 1);
    ind_last = length(summed_values) - find(flip(summed_values) > 0, 1) + 1;

    % trim data
    if numDim(data) == 1
        data = data(ind_first:ind_last);
    else
        switch dim_index
            case 1
                data = data(ind_first:ind_last, :, :);
                ind(1:2) = [ind_first, ind_last];
            case 2
                data = data(:, ind_first:ind_last, :);
                ind(3:4) = [ind_first, ind_last];
            case 3
                data = data(:, :, ind_first:ind_last);
                ind(5:6) = [ind_first, ind_last];
        end
    end
    
end