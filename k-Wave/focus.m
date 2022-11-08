function signal_mat = focus(kgrid, input_signal, source_mask, focus_position, sound_speed)
%FOCUS Create input signal based on source mask and focus position.
%
% DESCRIPTION:
%     focus takes a single input signal and a source mask and creates an
%     input signal matrix (with one input signal for each source point).
%     The appropriate time delays required to focus the signals at a given
%     position in Cartesian space are automatically added based on the user
%     inputs for focus_position and sound_speed.
%
% USAGE:
%     input_signal_mat = focus(kgrid, input_signal, source_mask, focus_position, sound_speed)
%
% INPUTS:
%     kgrid            - k-Wave grid object returned by kWaveGrid
%     input_signal     - single time series input
%     source_mask      - matrix specifying the positions of the time
%                        varying source distribution (i.e., source.p_mask
%                        or source.u_mask)  
%     focus_position   - position of the focus in Cartesian coordinates
%     sound_speed      - scalar sound speed
%
% OUTPUTS:
%     input_signal_mat - matrix of time series following the source points
%                        using MATLAB's column-wise linear index ordering  
%       
% ABOUT:
%     author           - Bradley Treeby
%     date             - 21st February 2012
%     last update      - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

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

% check that kgrid.t_array is defined
if strcmp(kgrid.t_array, 'auto')
    error('kgrid.t_array must be defined.');
end

% check that the input sound speed is scalar
if numel(sound_speed) ~= 1
    error('sound_speed input must be scalar.');
end

% calculate the distance from every point in the source mask to the focus
% position
switch kgrid.dim
    case 1
        dist = abs(kgrid.x(source_mask == 1) - focus_position(1));
    case 2
        dist = sqrt( (kgrid.x(source_mask == 1) - focus_position(1)).^2 + ...
            (kgrid.y(source_mask == 1) - focus_position(2)).^2 );
    case 3
        dist = sqrt( (kgrid.x(source_mask == 1) - focus_position(1)).^2 + ...
            (kgrid.y(source_mask == 1) - focus_position(2)).^2 + ...
            (kgrid.z(source_mask == 1) - focus_position(3)).^2 );
end

% convert the distance to units of time points
dist = round(dist ./ (kgrid.dt * sound_speed));

% convert time points to relative delays
dist = -(dist - max(dist(:)));
max_delay = max(dist(:));

% create an input matrix
signal_mat = zeros(length(dist), length(input_signal) + max_delay);

% assign the input signal
for source_index = 1:length(dist)
    delay = dist(source_index);
    signal_mat(source_index, :) = [zeros(1, delay), input_signal, zeros(1, max_delay - delay)];
end