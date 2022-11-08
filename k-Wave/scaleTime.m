function time = scaleTime(seconds)
%SCALETIME Convert seconds to hours, minutes, and seconds.
%
% DESCRIPTION:
%     scaleTime converts an integer number of seconds into hours, minutes,
%     and seconds, and returns a string with this information. 
%
% USAGE:
%     time = scaleTime(seconds)
%
% INPUTS:
%     seconds         - number of seconds
%
% OUTPUTS:
%     time            - string of scaled time   
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 16th July 2009
%     last update     - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also scaleSI

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

% switch to calculating years, weeks, and days if larger than 100 hours
if seconds > (60 * 60 * 100)
    years = floor(seconds / (60 * 60 * 24 * 365));
    seconds = seconds - years * 60 * 60 * 24 * 365;
    weeks = floor(seconds / (60 * 60 * 24 * 7));
    seconds = seconds - weeks * 60 * 60 * 24 * 7;
    days = floor(seconds / (60 * 60 * 24));
    seconds = seconds - days * 60 * 60 * 24;
else
    years = 0;
    weeks = 0;
    days = 0;
end
    
% calculate hours, minutes, and seconds
hours = floor(seconds / (60 * 60));
seconds = seconds - hours * 60 * 60;
minutes = floor( seconds / 60 );
seconds = seconds - minutes * 60;

% write out as a string, to keep the output manageable, only the largest
% three units are written
if years > 0
    time = [num2str(years) ' years, ' num2str(weeks) ' weeks, and ' num2str(days) ' days'];
elseif weeks > 0
    time = [num2str(weeks) ' weeks, ' num2str(days) ' days, and ' num2str(hours) ' hours'];
elseif days > 0
    time = [num2str(days) ' days, ' num2str(hours) ' hours, and ' num2str(minutes) ' min'];
elseif hours > 0
    time = [num2str(hours) 'hours ' num2str(minutes) 'min ' num2str(seconds) 's'];
elseif minutes > 0
    time = [num2str(minutes) 'min ' num2str(seconds) 's']; 
else
    time = [num2str(seconds) 's']; 
end