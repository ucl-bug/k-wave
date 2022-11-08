function path = getkWavePath(folder_name)
%GETKWAVEPATH Return pathname to the k-Wave Toolbox.
%
% DESCRIPTION:
%     getkWavePath returns the full directory pathname to the root 
%     directory in the k-Wave Toolbox using the slash direction native to
%     the users operating system. If a folder_name string is provided, this
%     is appended to the pathname.
%
% USAGE:
%     path = getkWavePath()
%     path = getkWavePath(folder_name)
%
% OPTIONAL INPUTS:
%     folder_name - folder name string to append to the pathname
%
%     Note: folder_name is not checked for existance, the string is
%     simply appended to the pathname with a trailing slash.
%
% OUTPUTS:
%     path        - full pathname to the root directory in the k-Wave
%                   Toolbox 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 30th June 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

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

% get the full pathname of the toolbox directory
path = mfilename('fullpath');

% trim off the m-file name
path = path(1:end - 12);

% add a folder name if given
if nargin ~= 0
    
    % extract required slash
    slash = path(end);
    
    % add folder name
    path = [path, folder_name, slash];
    
end