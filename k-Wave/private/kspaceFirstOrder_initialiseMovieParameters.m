% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to
%     initialise movie parameters. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 21st February 2011
%     last update - 19th February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% force getframe compatability with dual monitors
movegui(img);

% create new VideoWriter object (this is supported from MATLAB 2010b)
video_obj = VideoWriter(movie_name, movie_profile);

% adjust settings if specified by the user
if ~isempty(movie_args)
    for input_index = 1:2:length(movie_args)
        eval(['video_obj.' movie_args{input_index} ' = movie_args{input_index + 1};']);
    end
end
        
% open the object
open(video_obj);