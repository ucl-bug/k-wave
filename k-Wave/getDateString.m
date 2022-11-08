function date_string = getDateString()
%GETDATESTRING Create a string of the current date and time.
%
% DESCRIPTION:
%     getDateString returns a string of the current date and time using
%     datestr but replacing the white space and : character with an en
%     dash.
%
% USAGE:
%     date_string = getDateString()
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 14th October 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also clock, date, datestr

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

% get the current time
date_string = datestr(clock);

% replace the space and : characters with -
date_string = strrep(date_string, ' ', '-');
date_string = strrep(date_string, ':', '-');