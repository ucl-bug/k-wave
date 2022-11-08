function checkFieldNames(structure, field_names)
%CHECKFIELDNAMES Checks structure field names for existance.
%
% DESCRIPTION:
%     checkFieldNames checks the field names of a MATLAB structure against
%     a cell array of allowable field names. If a field name is found
%     within the structure that is not contained within the cell array, an
%     error is thrown with the field name given in the error message. The
%     function is useful for checking for spelling errors in field names
%     for structure inputs.
%
% USAGE:
%     checkFieldNames(structure, field_names)
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

% get list of field names
names = fieldnames(structure);

% loop through list
for names_index = 1:length(names)
    
    % check if current field name is in the list of allowable field names
    name_found = false;
    for field_index = 1:length(field_names)
        if strcmp(names{names_index}, field_names{field_index})
            name_found = true;
        end
    end
    
    % if name isn't found, throw error
    if ~name_found
    	error([names{names_index} ' is not a valid field for the structure ' inputname(1) '.']);
    end
    
end