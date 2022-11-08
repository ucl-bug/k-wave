function writeGrid(filename, grid_size, grid_spacing, pml_size, pml_alpha, Nt, dt, c_ref)  %#ok<*INUSD>
%WRITEGRID Write grid and PML properties to a k-Wave HDF5 file.
%
% DESCRIPTION:
%     writeGrid creates and writes the wavenumber grids and PML variables
%     required by the k-Wave C++ code to the HDF5 file specified by the
%     user. 
%
%     List of parameters that are written:
%         Nx
%         Ny
%         Nz
%         Nt
%         dt
%         dx
%         dy
%         dz
%         c_ref
%         pml_x_alpha
%         pml_y_alpha
%         pml_z_alpha
%         pml_x_size
%         pml_y_size
%         pml_z_size
%
% USAGE:
%     writeGrid(filename, grid_size, grid_spacing, pml_size, pml_alpha, Nt, dt, c_ref) 
%
% INPUTS:
%     filename            - filename and location of the input HDF5 file
%     grid_size           - [Nx, Ny, Nz]
%     grid_spacing        - [dx, dy, dz]
%     pml_size            - [pml_x_size, pml_y_size, pml_z_size]
%     pml_alpha           - [pml_x_alpha, pml_y_alpha, pml_z_alpha]
%     Nt                  - number of time points
%     dt                  - time step
%     c_ref               - scalar sound speed used in the k-space operator
%                           and to define the pml variables 
%
% ABOUT:
%     author              - Bradley Treeby
%     date                - 30th May 2013
%     last update         - 26th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2019 Bradley Treeby
%
% See also h5writeatt, writeAttributes, writeFlags, writeMatrix

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

%#ok<*INUSL>
%#ok<*NASGU>

% get literals
getH5Literals;

% unpack grid size inputs to make code easier to read
Nx          = grid_size(1);
Ny          = grid_size(2);
Nz          = grid_size(3);
dx          = grid_spacing(1);
dy          = grid_spacing(2);
dz          = grid_spacing(3);
pml_x_size  = pml_size(1);
pml_y_size  = pml_size(2);
pml_z_size  = pml_size(3);
pml_x_alpha = pml_alpha(1);
pml_y_alpha = pml_alpha(2);
pml_z_alpha = pml_alpha(3);

% =========================================================================
% STORE FLOATS
% =========================================================================

% list of variables stored as floats
variable_names = {...
    'dt', 'dx', 'dy', 'dz', ...
    'pml_x_alpha', 'pml_y_alpha', 'pml_z_alpha', ...
    'c_ref'};

% change float variables to be in single precision (float in C++), then
% add to HDF5 file
for index = 1:length(variable_names)

    % cast matrix to single precision
    eval([variable_names{index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% =========================================================================
% STORE INTEGERS
% =========================================================================

% integer variables
variable_names = {'Nx', 'Ny', 'Nz', 'Nt',...
    'pml_x_size' , 'pml_y_size' , 'pml_z_size'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end