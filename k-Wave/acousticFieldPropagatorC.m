function varargout = acousticFieldPropagatorC(amp_in, phase_in, dx, f0, c0, varargin)
%ACOUSTICFIELDPROPAGATOR Calculate acoustic field for CW source using C++ code.
%
% DESCRIPTION:
%     acousticFieldPropagatorC provides an interface to the C++ version of
%     acousticFieldPropagator (called acousticFieldPropagator-OMP). The
%     function works by appending the optional input 'SaveToDisk' to the
%     user inputs and then calling acousticFieldPropagator to save the
%     input files to disk. The C++ code is run using the system command.
%     The output files are then automatically loaded from disk and returned
%     in the same fashion as acousticFieldPropagator. The input and output
%     files are saved to the temporary directory native to the operating
%     system, and are deleted after the function runs.  
%
%     This function requires the C++ binary/executable of 
%     acousticFieldPropagator-OMP to be downloaded from 
%     http://www.k-wave.org/download.php and placed in the "binaries"
%     directory of the k-Wave toolbox. Alternatively, the name and
%     location of the binary can be specified using the optional input
%     parameters 'BinaryName' and 'BinariesPath'.
%
% USAGE:
%     see acousticFieldPropagator. 
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings.
%
%     'BinaryName'    - Name of the binary file (default =
%                       acousticFieldPropagator-OMP on linux, 
%                       acousticFieldPropagator-OMP.exe on windows).
%     'BinaryPath'    - Path of the binary file (default = binaries/).
%     'DataName'      - Prefix used to generate a custom name for the
%                       input and output data files (this is appended
%                       with _input.h5 and _output.h5) (default =
%                       kwave_<input/output>_data_<date>.h5).
%     'DataPath'      - Location of the folder where the input and output
%                       HDF5 files should be stored (default = tempdir).
%     'DeleteData'    - Boolean controlling whether the input and output
%                       HDF5 files should be deleted after running the
%                       simulation (default = true).
%
% ABOUT:
%     author          - Bradley Treeby and Jakub Budisky
%     date            - 27th February 2017
%     last update     - 4th April 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Bradley Treeby and Jakub Budisky
%
% See also acousticFieldPropagator

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

%#ok<*AGROW>

% start time
global_start_time = clock;

% check the inputs are 3D
if numDim(amp_in) ~= 3
    error('Only 3D simulations are supported using acousticFieldPropagatorC.');
end

% set options string to return output pressure in complex (interleaved)
% form
options_string = ' -c ';

% extract the optional input arguments
if nargin > 5
    input_args = varargin;
else
    input_args = {};
end

% check for a user defined location for the binary
if any(strcmp('BinaryPath', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('BinaryPath', input_args));
    
    % read the value of the parameter
    binary_path = input_args{index + 1};
    
    % check for a trailing slash
    if ~strcmp(binary_path(end), filesep)
        binary_path = [binary_path filesep];
    end
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default location to the path of the m-file
    binary_path = mfilename('fullpath');
    
    % remove name of m-file, and append binaries folder
    binary_path = [binary_path(1:end - length(mfilename)), 'binaries', filesep];
    
end
    
% check for a user defined name for the binary
if any(strcmp('BinaryName', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('BinaryName', input_args));
    
    % read the value of the parameter
    binary_name = input_args{index + 1};
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default name for the binary
    if isunix
        binary_name = 'acousticFieldPropagator-OMP';
    else
        binary_name = 'acousticFieldPropagator-OMP.exe';
    end
    
end

% check the binaries exist and are in the correct place before doing
% anything else
if ~exist([binary_path, binary_name], 'file')
    
    error(['The binary file ' binary_name ' could not be found in ' ...
        binary_path '. To use the C++ code, the C++ binaries for your ' ...
        'operating system must be downloaded from ' ...
        'www.k-wave.org/download.php and placed in the binaries folder.']);
    
end

% check for a user defined location for the input and output files
if any(strcmp('DataPath', input_args))
   
    % find the location of the optional input in the list
    index = find(strcmp('DataPath', input_args));
    
    % read the value of the parameter
    data_path = input_args{index + 1};
    
    % check for a trailing slash
    if ~strcmp(data_path(end), filesep)
        data_path = [data_path filesep];
    end
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default path
    data_path = tempdir;
    
end
    
% check for a user defined name for the input and output files
if any(strcmp('DataName', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('DataName', input_args));
    
    % read the value of the parameter
    name_prefix = input_args{index + 1};
    input_filename = [name_prefix '_input.h5'];
    output_filename = [name_prefix '_output.h5'];
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];    
    
else
    
    % set the filename inputs to store data in the default temp directory
    date_string = getDateString;
    input_filename = ['kwave_input_data' date_string '.h5'];
    output_filename = ['kwave_output_data' date_string '.h5'];
    
end

% add pathname to input and output filenames
input_filename = [data_path, input_filename];
output_filename = [data_path, output_filename];

% check for delete data input
if any(strcmp('DeleteData', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('DeleteData', input_args));
    
    % read the value of the parameter
    delete_data = input_args{index + 1};
    
    % check option is true or false
    if ~islogical(delete_data)
        error('Optional input ''DeleteData'' must be Boolean');
    end
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set data to be deleted
    delete_data = true;
    
end

% run the MATLAB code first to generate the input file and save to disk
acousticFieldPropagator(amp_in, phase_in, dx, f0, c0, 'SaveToDisk', input_filename, input_args{:});

% store time to save input
compute_time.save_input = etime(clock, global_start_time);

% start timer
compute_start_time = clock;

% run the simulation in C++ and print outputs to the MATLAB command line
if isunix
    
    % prepend spaces in linux pathnames with \ to allow cd to work
    binary_path = strrep(binary_path, ' ', '\ ');
    
    % clear the library path to prevent a conflict with the FFTW libraries
    % loaded automatically by MATLAB, and run linux binary
    run_string = ['system(''export LD_LIBRARY_PATH=; cd ' binary_path '; ./' binary_name ' -i ' input_filename ' -o ' output_filename ' ' options_string ''' ,''-echo'');'];
    eval(run_string);
    
else
    
    % run Windows binary, placing the input and output filenames in double
    % quotations to avoid problems with spaces
    run_string = ['system(''cd /d ' binary_path ' & ' binary_name ' -i "' input_filename '" -o "' output_filename '" ' options_string ''' ,''-echo'');'];
    eval(run_string);
    
end

% store time to run simulation
compute_time.simulation = etime(clock, compute_start_time);

% start timer
load_start_time = clock;

% load the C++ data back from disk using h5read
pressure = h5read(output_filename, '/pressure_out');

% reformat the data to MATLAB complex format
pressure = pressure(:, :, :, 1) + 1i .* pressure(:, :, :, 2);

% compute amplitude and phase from the two beam patterns
if nargout > 1
    amp_out   = abs(pressure);
    phase_out = angle(pressure);
end

% delete the input and output files
if delete_data
    delete(input_filename);
    delete(output_filename);
end

% store time to load output
compute_time.load_output = etime(clock, load_start_time);

% store total time
compute_time.total = etime(clock, global_start_time);

% assign outputs
switch nargout
    case 1
        varargout{1} = pressure;
    case 2
        varargout{1} = amp_out;
        varargout{2} = phase_out;
    case 3
        varargout{1} = pressure;
        varargout{2} = amp_out;
        varargout{3} = phase_out;
    case 4
        varargout{1} = pressure;
        varargout{2} = amp_out;
        varargout{3} = phase_out;
        varargout{4} = compute_time;
end

% Here, compute_time is an undocumented output option to return a timing
% structure containing the following fields:
%     .total
%     .save_input
%     .simulation
%     .load_output