function kgrid = makeGrid(varargin)
%MAKEGRID Create k-Wave grid structure.
%
% DESCRIPTION:
%     See kWaveGrid for inputs. Note, this function will be deprecated in
%     a future version of k-Wave.
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 12th March 2009
%     last update   - 8th June 2017
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

% display warning
disp('WARNING: makeGrid will be deprecated in a future version of k-Wave.')
disp('         Update codes to use the syntax kgrid = kWaveGrid(...).');

try
    
    % create a new instance of the kWaveGrid class
    kgrid = kWaveGrid(varargin{:});
    
catch %#ok<CTCH>
    
    % throw error if user defined classes are not supported
    error('User defined classes are not supported in this version of MATLAB. Please use a more recent version of MATLAB, or k-Wave V1.1.');
    
end