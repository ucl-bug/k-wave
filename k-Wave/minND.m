function [min_val, ind] = minND(mat)
%MAXND Return the value and indices of the smallest value in an N-D array.
%
% DESCRIPTION:
%     minND returns the smallest value in the N-D array given by mat. The
%     indices of the smallest value within the array can also be
%     returned, where ind = [ind_dim_1, ind_dim_2, ... , ind_dim_N]. The
%     indices can be used to directly address the smallest value, i.e.,
%     min_val = mat(ind).
%
% USAGE:
%     min_val = minND(mat)
%     [min_val, ind] = minND(mat)
%
% INPUTS:
%     mat             - input array
%
% OUTPUTS:
%     min_val         - smallest value in mat
%     ind             - array of indices corresponding to the location of
%                       the smallest value in mat
%       
% ABOUT:
%     author          - Bradley Treeby
%     date            - 3rd April 2012
%     last update     - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

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

% get maximum value
[min_val, linear_index] = min(mat(:));
    
% convert linear index to matrix indices
if nargout == 2
    [ind{1:ndims(mat)}] = ind2sub(size(mat), linear_index);
    ind = cell2mat(ind);
end