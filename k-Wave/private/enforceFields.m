function enforceFields(structure, field_names)
%ENFORCEFIELDS Check structure field names for existance.
%
% DESCRIPTION:
%     enforceFields checks a MATLAB structure for the existance of a set of
%     required field name defined by a cell array. If a field name within
%     the cell array is not found within the structure, an error is thrown
%     with the name of the missing field given in the error message. The
%     function is useful for enforcing a set of required fields within a
%     structure. 
%
% USAGE:
%     enforceFields(structure, field_names)
%
% INPUTS:
%     structure   - MATLAB structure
%     field_names - cell array of allowable field names
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 13th October 2009
%     last update - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also fieldnames

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

% get list of field names in the structure
names = fieldnames(structure);

% loop through list of require field names
for field_index = 1:length(field_names)
    
    % check if each required field names is in the structure
    field_found = false;
    for names_index = 1:length(names)
        if strcmp(field_names{field_index}, names{names_index})
            field_found = true;
        end
    end
    
    % if the field is missing, throw an error
    if ~field_found
    	error(['The field ' inputname(1) '.' field_names{field_index} ' must be defined.']);
    end
    
end