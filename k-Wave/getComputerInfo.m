function comp_info = getComputerInfo()
%GETCOMPUTERINFO Return information about computer and k-Wave version.
%
% DESCRIPTION:
%     getComputerInfo returns information about the computer currently
%     being used, including the versions of MATLAB and k-Wave.
%
% USAGE:
%     comp_info = getComputerInfo()
%
% OUTPUTS:
%     comp_info   - MATLAB structure containing the following fields:
%                       .date
%                       .computer_name
%                       .operating_system_type
%                       .operating_system
%                       .user_name
%                       .matlab_version
%                       .kwave_version
%                       .kwave_path
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 19 February 2014
%     last update - 21 April 2021
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2021 Bradley Treeby
%
% See also ver

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

% get date
comp_info.date = date;

% get computer name, removing carriage returns
[~, hname] = system('hostname');
comp_info.computer_name = regexprep(hname,'\r\n|\n|\r','');

% get os
comp_info.operating_system_type = computer;
comp_info.operating_system = system_dependent('getos');

% get username
if isunix() 
    comp_info.user_name = getenv('USER'); 
else 
    comp_info.user_name = getenv('username'); 
end

% get matlab version
v = ver('matlab');
comp_info.matlab_version = [v.Version ' ' v.Release];

% get kwave version
comp_info.kwave_version = getkWaveVersion;

% get kwave path
comp_info.kwave_path = getkWavePath;