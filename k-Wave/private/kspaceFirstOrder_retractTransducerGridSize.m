% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to retract
%     the grid size in an object of the kWaveTransducer class after
%     simulation with 'PMLInside', false. This is needed as kWaveTransducer
%     is implemented as a handle class.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 26th September 2012
%     last update - 18th December 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2018 Bradley Treeby

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

% resize the transducer object if the grid has been expanded
if ~flags.pml_inside && (isa(source, 'kWaveTransducer') || isa(sensor, 'kWaveTransducer'))
    
    % retract by the pml size if retract size not defined
    if ~exist('retract_size', 'var')
        retract_size = [pml_x_size, pml_y_size, pml_z_size]; 
    end
    
    % check if the sensor is a transducer
    if isa(sensor, 'kWaveTransducer')
        
        % retract the transducer mask
        sensor.retract_grid(retract_size);
        
    end
        
    % check if the source is a transducer, and if so, and different
    % transducer to the sensor 
    if isa(source, 'kWaveTransducer') && ~(isa(sensor, 'kWaveTransducer') && isequal(sensor, source))
        
        % retract the transducer mask
        source.retract_grid(retract_size);
        
    end
    
end