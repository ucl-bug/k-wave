function writeFlags(filename) 
%WRITEGRID Write input flags to a k-Wave HDF5 file.
%
% DESCRIPTION:
%     writeFlags reads the input HDF5 file and derives and writes the
%     required source and medium flags based on the datasets present in the
%     file. For example, if the file contains a data set named 'BonA', the
%     nonlinear_flag will be written as true. Conditional flags are also
%     written. The source mode flags are written when appropriate if they
%     are not already present in the file. The default source mode is
%     'additive'. 
%
%     List of flags that are always written
%         ux_source_flag
%         uy_source_flag
%         uz_source_flag
%         sxx_source_flag
%         sxy_source_flag
%         sxz_source_flag
%         syy_source_flag
%         syz_source_flag
%         szz_source_flag
%         p_source_flag
%         p0_source_flag
%         transducer_source_flag
%         nonuniform_grid_flag
%         nonlinear_flag
%         absorbing_flag
%         axisymmetric_flag
%         elastic_flag
%         sensor_mask_type
%
%     List of conditional flags
%         u_source_mode
%         u_source_many
%         p_source_mode 
%         p_source_many
%         s_source_mode 
%         s_source_many
%
% USAGE:
%     writeFlags(filename) 
%
% INPUTS:
%     filename            - filename and location of the input HDF5 file
%
% ABOUT:
%     author              - Bradley Treeby
%     date                - 31st May 2013
%     last update         - 26th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2018 Bradley Treeby
%
% See also h5writeatt, writeAttributes, writeGrid, writeMatrix

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

%#ok<*NASGU>

% get literals
getH5Literals;

% list of variable names
variable_names = {...
    'ux_source_flag', 'uy_source_flag', 'uz_source_flag', ...
    'sxx_source_flag', 'syy_source_flag', 'szz_source_flag', ...
    'sxy_source_flag', 'sxz_source_flag', 'syz_source_flag', ...
    'p_source_flag', ...
    'p0_source_flag', ...
    'transducer_source_flag', ...
    'nonuniform_grid_flag', ...
    'nonlinear_flag', ...
    'absorbing_flag', ...
    'axisymmetric_flag', ...
    'elastic_flag', ...
    'sensor_mask_type'};

% load metadata from HDF5 file
filename_info = h5info(filename);

% extract the names of the stored datasets
names = {filename_info.Datasets.Name};

% --------------------
% u source
% --------------------

% check for ux_source_input and set ux_source_flag
if any(strcmp(names, 'ux_source_input'))
    
    % find index of ux_source_input in names and readout size
    index = find(strcmp(names, 'ux_source_input') == 1);
    ux_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        u_source_many = 0;
    else
        u_source_many = 1;
    end
    
else
    ux_source_flag = 0;
end

% check for uy_source_input and set uy_source_flag
if any(strcmp(names, 'uy_source_input'))
    
    % find index of uy_source_input in names and readout size
    index = find(strcmp(names, 'uy_source_input') == 1);
    uy_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        u_source_many = 0;
    else
        u_source_many = 1;
    end
    
else
    uy_source_flag = 0;
end

% check for uz_source_input and set uz_source_flag
if any(strcmp(names, 'uz_source_input'))
    
    % find index of uz_source_input in names and readout size
    index = find(strcmp(names, 'uz_source_input') == 1);
    uz_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        u_source_many = 0;
    else
        u_source_many = 1;
    end
    
else
    uz_source_flag = 0;
end

% add conditional flag to list of variable names
if exist('u_source_many') %#ok<EXIST>
    variable_names = [variable_names, {'u_source_many'}];
end

% write u_source mode if not already in file (1 is Additive, 0 is Dirichlet)
if (ux_source_flag || uy_source_flag || uz_source_flag) && ~any(strcmp(names, 'u_source_mode'))
    u_source_mode = 1;
    variable_names = [variable_names, {'u_source_mode'}];
end

% --------------------
% s source
% --------------------

% check for sxx_source_input and set sxx_source_flag
if any(strcmp(names, 'sxx_source_input'))
    
    % find index of sxx_source_input in names and readout size
    index = find(strcmp(names, 'sxx_source_input') == 1);
    sxx_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    sxx_source_flag = 0;
end

% check for syy_source_input and set syy_source_flag
if any(strcmp(names, 'syy_source_input'))
    
    % find index of syy_source_input in names and readout size
    index = find(strcmp(names, 'syy_source_input') == 1);
    syy_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    syy_source_flag = 0;
end

% check for szz_source_input and set szz_source_flag
if any(strcmp(names, 'szz_source_input'))
    
    % find index of szz_source_input in names and readout size
    index = find(strcmp(names, 'szz_source_input') == 1);
    szz_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    szz_source_flag = 0;
end

% check for sxy_source_input and set sxy_source_flag
if any(strcmp(names, 'sxy_source_input'))
    
    % find index of sxy_source_input in names and readout size
    index = find(strcmp(names, 'sxy_source_input') == 1);
    sxy_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    sxy_source_flag = 0;
end

% check for sxz_source_input and set sxz_source_flag
if any(strcmp(names, 'sxz_source_input'))
    
    % find index of sxz_source_input in names and readout size
    index = find(strcmp(names, 'sxz_source_input') == 1);
    sxz_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    sxz_source_flag = 0;
end

% check for syz_source_input and set syz_source_flag
if any(strcmp(names, 'syz_source_input'))
    
    % find index of syz_source_input in names and readout size
    index = find(strcmp(names, 'syz_source_input') == 1);
    syz_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        s_source_many = 0;
    else
        s_source_many = 1;
    end
    
else
    syz_source_flag = 0;
end

% add conditional flag to list of variable names
if exist('s_source_many') %#ok<EXIST>
    variable_names = [variable_names, {'s_source_many'}];
end

% write s_source mode if not already in file (1 is Additive, 0 is Dirichlet)
if (sxx_source_flag || syy_source_flag || szz_source_flag || sxy_source_flag || sxz_source_flag || syz_source_flag) && ~any(strcmp(names, 's_source_mode'))
    s_source_mode = 1;
    variable_names = [variable_names, {'s_source_mode'}];
end

% --------------------
% p source
% --------------------

% check for p_source_input and set p_source_flag
if any(strcmp(names, 'p_source_input'))
    
    % find index of p_source_input in names and readout size
    index = find(strcmp(names, 'p_source_input') == 1);
    p_source_flag = filename_info.Datasets(index).Dataspace.Size(2);
    
    % check for multiple sources
    if filename_info.Datasets(index).Dataspace.Size(1) == 1
        p_source_many = 0;
    else
        p_source_many = 1; 
    end
    
    % add conditional flag to list of variable names
    variable_names = [variable_names, {'p_source_many'}];
    
else
    p_source_flag = 0;
end

% write p_source mode if not already in file (1 is Additive, 0 is Dirichlet)
if p_source_flag && ~any(strcmp(names, 'p_source_mode'))
    p_source_mode = 1;
    variable_names = [variable_names, {'p_source_mode'}];
end

% check for p0_source_input and set p0_source_flag
p0_source_flag = any(strcmp(names, 'p0_source_input'));
    
% --------------------
% additional flags
% --------------------

% check for transducer_source_input and set transducer_source_flag
transducer_source_flag = any(strcmp(names, 'transducer_source_input'));

% check for BonA and set nonlinear flag
nonlinear_flag = any(strcmp(names, 'BonA'));

% check for alpha_coeff and set absorbing flag
absorbing_flag = any(strcmp(names, 'alpha_coeff'));

% check for lambda and set elastic flag
elastic_flag = any(strcmp(names, 'lambda'));

% set axisymmetric grid flag to false
axisymmetric_flag = 0;

% set nonuniform grid flag to false
nonuniform_grid_flag = 0;

% check for sensor_mask_index and sensor_mask_corners
if any(strcmp(names, 'sensor_mask_index'))
    sensor_mask_type = 0;
elseif any(strcmp(names, 'sensor_mask_corners'))
    sensor_mask_type = 1;
else
    error('Either sensor_mask_index or sensor_mask_corners must be defined in the input file');
end

% --------------------
% write flags to file
% --------------------

% change all the index variables to be in 64-bit unsigned integers (long in
% C++) and write to file
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end