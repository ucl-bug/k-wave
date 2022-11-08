function writeAttributes(filename, file_description)
%WRITEMATRIX Write attributes to a k-Wave HDF5 file.
%
% DESCRIPTION:
%     writeAttributes writes a set of additional attributes to the HDF5
%     file specified by the user. These attributes specify the file
%     version, k-Wave version, type, and so on, and are required by the
%     k-Wave C++ simulation code.
%
%     If the file_description is not specified, a description is generated
%     based on the computer username and MATLAB version.
%
% USAGE:
%     writeAttributes(filename)
%     writeAttributes(filename, file_description)
%
% INPUTS:
%     filename            - name of HDF5 file to write matrix to
%     file_description    - custom file description
%
% ABOUT:
%     author              - Bradley Treeby and Jiri Jaros
%     date                - 20th February 2013
%     last update         - 21st October 2022
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2022 Bradley Treeby and Jiri Jaros
%
% See also h5writeatt, writeFlags, writeGrid, writeMatrix

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

% get literals
getH5Literals;

% get computer infor
comp_info = getComputerInfo;

% set file description if not provided by user
if nargin == 1
    file_description = ['Input data created by ' comp_info.user_name ' running MATLAB ' comp_info.matlab_version ' on ' comp_info.operating_system_type];
end

% set additional file attributes
if verLessThan('matlab', '9.8')
    h5writeatt(filename, '/', FILE_MAJOR_VER_ATT_NAME,     HDF_FILE_MAJOR_VERSION);
    h5writeatt(filename, '/', FILE_MINOR_VER_ATT_NAME,     HDF_FILE_MINOR_VERSION);
    h5writeatt(filename, '/', CREATED_BY_ATT_NAME,         ['k-Wave ' getkWaveVersion]);
    h5writeatt(filename, '/', FILE_DESCR_ATT_NAME,         file_description);
    h5writeatt(filename, '/', FILE_TYPE_ATT_NAME,          HDF_INPUT_FILE);
    h5writeatt(filename, '/', FILE_CREATION_DATE_ATT_NAME, getDateString);
else
    h5writeatt(filename, '/', FILE_MAJOR_VER_ATT_NAME,     HDF_FILE_MAJOR_VERSION,      'TextEncoding', 'system');
    h5writeatt(filename, '/', FILE_MINOR_VER_ATT_NAME,     HDF_FILE_MINOR_VERSION,      'TextEncoding', 'system');
    h5writeatt(filename, '/', CREATED_BY_ATT_NAME,         ['k-Wave ' getkWaveVersion], 'TextEncoding', 'system');
    h5writeatt(filename, '/', FILE_DESCR_ATT_NAME,         file_description,            'TextEncoding', 'system');
    h5writeatt(filename, '/', FILE_TYPE_ATT_NAME,          HDF_INPUT_FILE,              'TextEncoding', 'system');
    h5writeatt(filename, '/', FILE_CREATION_DATE_ATT_NAME, getDateString,               'TextEncoding', 'system');
end