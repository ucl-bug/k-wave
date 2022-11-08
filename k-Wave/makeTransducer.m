function transducer = makeTransducer(kgrid, transducer_properties)
%MAKETRANSDUCER Create k-Wave ultrasound transducer.
%
% DESCRIPTION:
%       See kWaveTransducer for inputs. Note, this function will be
%       deprecated in a future version of k-Wave.
%
% ABOUT:
%       author          - Bradley Treeby
%       date            - 28th July 2011
%       last update     - 4th June 2017
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

% display warning
disp('WARNING: makeTransducer will be deprecated in a future version of k-Wave.')
disp('         Update codes to use the syntax transducer = kWaveTransducer(...).');

try
    
    % create a new instance of the kWaveTransducer class
    transducer = kWaveTransducer(kgrid, transducer_properties);
    
catch ME
    
    % if user defined classes aren't supported, throw an error 
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        error('The transducer cannot be created because user defined classes are not supported in your version of MATLAB. To use this functionality, please try using a newer MATLAB version.');
    end
    rethrow(ME);
    
end
