function mat_new = expandMatrix(mat, exp_coeff, edge_val)
%EXPANDMATRIX Enlarge a matrix by extending the edge values.
%
% DESCRIPTION:
%     expandMatrix enlarges an input matrix by extension of the values at
%     the outer faces of the matrix (endpoints in 1D, outer edges in 2D,
%     outer surfaces in 3D). Alternatively, if an input for edge_val is
%     given, all expanded matrix elements will have this value. The values
%     for exp_coeff are forced to be real positive integers (or zero).
%     Note, indexing is done inline with other k-Wave functions using
%     mat(x) in 1D, mat(x, y) in 2D, and mat(x, y, z) in 3D.
%
%     For example, running 
%
%         mat = magic(3)
%         expandMatrix(mat, 1)
%         expandMatrix(mat, [2 0 1 0], 0)
%
%     will give the outputs 
%
%         ans =
% 
%             8     8     1     6     6
%             8     8     1     6     6
%             3     3     5     7     7
%             4     4     9     2     2
%             4     4     9     2     2
% 
%         ans =
% 
%             0     0     0     0
%             0     0     0     0
%             0     8     1     6
%             0     3     5     7
%             0     4     9     2
%
% USAGE:
%     mat_new = expandMatrix(mat, exp_coeff)
%     mat_new = expandMatrix(mat, exp_coeff, edge_val)
%
% INPUTS:
%     mat         - the matrix to enlarge
%     exp_coeff   - the number of elements to add in each dimension
%                   in 1D: [a] or [x_start, x_end]
%                   in 2D: [a] or [x, y] or 
%                          [x_start, x_end, y_start, y_end]
%                   in 3D: [a] or [x, y, z] or 
%                          [x_start, x_end, y_start, y_end, z_start, z_end]
%                          (here 'a' is applied to all dimensions)
%
% OPTIONAL INPUTS:
%     edge_val    - value to use in the matrix expansion
%
% OUTPUTS:
%     mat_new     - expanded matrix
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 13th October 2009
%     last update - 15 March 2022
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2022 Bradley Treeby
%
% See also resize, trimZeros

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

% check to see if a value for edge_val has been given
if nargin == 3
    extend_edges = false;
elseif nargin == 2
    extend_edges = true;
    edge_val = 1;
else
    error('Incorrect number of inputs.');
end

% force the values given for edge_val to be real positive integers
edge_val = round(abs(real(edge_val)));

% check the class of the input matrix
mat_class = class(mat);

% check if the input matrix is logical, if so, create the expanded matrix
% using int8 (1-byte per element) and then convert back to logical
if islogical(mat)
    mat = int8(mat);
    mat_class = 'int8';
    convert_to_logical = true;
else
    convert_to_logical = false;
end

switch numDim(mat)
    
    case 1
        
        % extract expansion coefficients
        if length(exp_coeff) == 2
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(2);
        elseif length(exp_coeff) == 1
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(1);
        else
            error('exp_coeff must be a 1 or 2 element array.');
        end
               
        % create a new expanded matrix
        x = length(mat);
        x_new = x + x1_size + x2_size;
        mat_new = edge_val*ones(x_new, 1, mat_class);
        
        % create indexes to allow the original matrix to be placed into the
        % expanded matrix
        x1 = x1_size + 1;
        x2 = x_new - x2_size;
        mat_new(x1:x2) = mat;
        
        % replace the remaining values by extending the edges
        if extend_edges
            mat_new(1:x1-1) = mat(1);
            mat_new(x2+1:end) = mat(end); 
        end
        
    case 2

        % extract expansion coefficients
        if length(exp_coeff) == 4
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(2);
            y1_size = exp_coeff(3);
            y2_size = exp_coeff(4);
        elseif length(exp_coeff) == 2
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(1);
            y1_size = exp_coeff(2);
            y2_size = exp_coeff(2);    
        elseif length(exp_coeff) == 1
            x1_size = exp_coeff;
            x2_size = exp_coeff;
            y1_size = exp_coeff;
            y2_size = exp_coeff;
        else
            error('exp_coeff must be a 1, 2, or 4 element array.');
        end

        % create a new expanded matrix
        [x, y] = size(mat);
        x_new = x + x1_size + x2_size;
        y_new = y + y1_size + y2_size;
        mat_new = edge_val*ones(x_new, y_new, mat_class);

        % create indexes to allow the original matrix to be placed into the
        % expanded matrix
        x1 = x1_size + 1;
        x2 = x_new - x2_size;
        y1 = y1_size + 1;
        y2 = y_new - y2_size;
        mat_new(x1:x2, y1:y2) = mat;        

        % replace the remaining values by extending the outer surfaces
        if extend_edges        
                        
            % extend the edge values
            mat_new(1:x1-1, y1:y2) = repmat(mat(1, :), x1_size, 1);
            mat_new(x2+1:end, y1:y2) = repmat(mat(end, :), x2_size, 1);
            mat_new(x1:x2, 1:y1-1) = repmat(mat(:, 1), 1, y1_size);
            mat_new(x1:x2, y2+1:end) = repmat(mat(:, end), 1, y2_size);
            
            % extend the corner values
            mat_new(1:x1-1, 1:y1-1) = mat(1, 1)*ones(x1_size, y1_size, mat_class);
            mat_new(1:x1-1, y2+1:end) = mat(1, end)*ones(x1_size, y2_size, mat_class);
            mat_new(x2+1:end, 1:y1-1) = mat(end, 1)*ones(x2_size, y1_size, mat_class);
            mat_new(x2+1:end, y2+1:end) = mat(end, end)*ones(x2_size, y2_size, mat_class);
        end    
                
    case 3
        
        % extract expansion coefficients
        if length(exp_coeff) == 6
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(2);
            y1_size = exp_coeff(3);
            y2_size = exp_coeff(4);            
            z1_size = exp_coeff(5);
            z2_size = exp_coeff(6);
        elseif length(exp_coeff) == 3
            x1_size = exp_coeff(1);
            x2_size = exp_coeff(1);
            y1_size = exp_coeff(2);
            y2_size = exp_coeff(2);              
            z1_size = exp_coeff(3);
            z2_size = exp_coeff(3);    
        elseif length(exp_coeff) == 1
            x1_size = exp_coeff;
            x2_size = exp_coeff;
            y1_size = exp_coeff;
            y2_size = exp_coeff;            
            z1_size = exp_coeff;
            z2_size = exp_coeff;
        else
            error('exp_coeff must be a 1, 3, or 6 element array.');
        end

        % create a new expanded matrix
        [x, y, z] = size(mat);
        x_new = x + x1_size + x2_size;
        y_new = y + y1_size + y2_size;
        z_new = z + z1_size + z2_size;
        mat_new = edge_val*ones(x_new, y_new, z_new, mat_class);

        % create indexes to allow the original matrix to be placed into the
        % expanded matrix
        x1 = x1_size + 1;
        x2 = x_new - x2_size;
        y1 = y1_size + 1;
        y2 = y_new - y2_size;        
        z1 = z1_size + 1;
        z2 = z_new - z2_size;
        mat_new(x1:x2, y1:y2, z1:z2) = mat; 

        % replace the remaining values by extending the outer surfaces
        if extend_edges        
                        
            % extend the face values
            mat_new(1:x1-1, y1:y2, z1:z2) = repmat(mat(1, :, :), [x1_size, 1, 1]);
            mat_new(x2+1:end, y1:y2, z1:z2) = repmat(mat(end, :, :), [x2_size, 1, 1]);
            mat_new(x1:x2, 1:y1-1, z1:z2) = repmat(mat(:, 1, :), [1, y1_size, 1]);
            mat_new(x1:x2, y2+1:end, z1:z2) = repmat(mat(:, end, :), [1, y2_size, 1]);
            mat_new(x1:x2, y1:y2, 1:z1-1) = repmat(mat(:, :, 1), [1, 1, z1_size]);
            mat_new(x1:x2, y1:y2, z2+1:end) = repmat(mat(:, :, end), [1, 1, z2_size]);            
            
            % extend the edge values
            mat_new(1:x1-1, 1:y1-1, z1:z2) = repmat(mat(1, 1, :), [x1_size, y1_size, 1]);
            mat_new(1:x1-1, y2+1:end, z1:z2) = repmat(mat(1, end, :), [x1_size, y2_size, 1]);
            mat_new(1:x1-1, y1:y2, 1:z1-1) = repmat(mat(1, :, 1), [x1_size, 1, z1_size]);
            mat_new(1:x1-1, y1:y2, z2+1:end) = repmat(mat(1, :, end), [x1_size, 1, z2_size]);
            mat_new(x2+1:end, 1:y1-1, z1:z2) = repmat(mat(end, 1, :), [x2_size, y1_size, 1]);
            mat_new(x2+1:end, y2+1:end, z1:z2) = repmat(mat(end, end, :), [x2_size, y2_size, 1]);
            mat_new(x2+1:end, y1:y2, 1:z1-1) = repmat(mat(end, :, 1), [x2_size, 1, z1_size]);
            mat_new(x2+1:end, y1:y2, z2+1:end) = repmat(mat(end, :, end), [x2_size, 1, z2_size]);
            mat_new(x1:x2, 1:y1-1, 1:z1-1) = repmat(mat(:, 1, 1), [1, y1_size, z1_size]);
            mat_new(x1:x2, y2+1:end, 1:z1-1) = repmat(mat(:, end, 1), [1, y2_size, z1_size]);
            mat_new(x1:x2, 1:y1-1, z2+1:end) = repmat(mat(:, 1, end), [1, y1_size, z2_size]);
            mat_new(x1:x2, y2+1:end, z2+1:end) = repmat(mat(:, end, end), [1, y2_size, z2_size]);
            
            % extend corner values
            mat_new(1:x1-1, 1:y1-1, 1:z1-1) = mat(1, 1, 1)*ones(x1_size, y1_size, z1_size, mat_class);
            mat_new(1:x1-1, y2+1:end, 1:z1-1) = mat(1, end, 1)*ones(x1_size, y2_size, z1_size, mat_class);
            mat_new(1:x1-1, 1:y1-1, z2+1:end) = mat(1, 1, end)*ones(x1_size, y1_size, z2_size, mat_class);
            mat_new(1:x1-1, y2+1:end, z2+1:end) = mat(1, end, end)*ones(x1_size, y2_size, z2_size, mat_class);
            mat_new(x2+1:end, 1:y1-1, 1:z1-1) = mat(end, 1, 1)*ones(x2_size, y1_size, z1_size, mat_class);
            mat_new(x2+1:end, y2+1:end, 1:z1-1) = mat(end, end, 1)*ones(x2_size, y2_size, z1_size, mat_class);
            mat_new(x2+1:end, 1:y1-1, z2+1:end) = mat(end, 1, end)*ones(x2_size, y1_size, z2_size, mat_class);
            mat_new(x2+1:end, y2+1:end, z2+1:end) = mat(end, end, end)*ones(x2_size, y2_size, z2_size, mat_class);

        end        
        
    otherwise
        error('Input matrix must be 1, 2 or 3 dimensional.');
end

if convert_to_logical
    mat_new = logical(mat_new);
end