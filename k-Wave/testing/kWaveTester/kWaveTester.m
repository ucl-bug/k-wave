function kWaveTester(options)
%KWAVETESTER Test all k-Wave options.
%
% DESCRIPTION:
%     kWaveTester runs a series of simulations using kspaceFirstOrder1D,
%     kspaceFirstOrder2D, and kspaceFirstOrder3D and kspaceFirstOrderAS to
%     test for runtime errors. For kspaceFirstOrder2D, kspaceFirstOrder3D,
%     and kspaceFirstOrderAS, the outputs can also be compared with the C++
%     and CUDA codes. The output log is saved to disk. 
%
% USAGE:
%     kWaveTester(options)
%
% INPUTS:
%
%     ---------------
%     GENERAL OPTIONS
%     ---------------
%
%     options.use_nonuniform_grid 
%         - Boolean to use non-uniform grid (NOT CURRENTLY USED)
%     options.force_plot_off 
%         - Boolean to force MATLAB k-Wave plotting to be off
%     options.output_folder
%         - path to folder to store logs and images
%     options.save_test_log
%         - Boolean to save the test log
%     options.split_log_after_n_tests
%         - integer after which to split the test log into a new file
%     options.continue_past_errors
%         - Boolean controlling whether kWaveTester should continue if an
%           error is encountered
%
%     ---------------
%     C++ OPTIONS
%     ---------------
%
%     options.run_cpp_comparison_tests 
%         - Boolean to run comparisons with C++ code
%     options.plot_cpp_errors
%         - Boolean to plot the errors against the C++ code
%     options.save_cpp_comparison_plots_to_disk
%         - Boolean to save the errors plots to disk
%     options.image_foldername
%         - path for folder to store error plots
%     options.use_gpu_code
%         - Boolean to call ...2DG or 3DG instead of 2DC or 3DC
%     options.cpp_binary_path
%         - path for C++ binary
%     options.cpp_binary_name
%         - name of the C++ binary to test
%     options.cpp_save_to_disk_only
%         - Boolean to just save the HDF5 input file and not run the test 
%           (only used if options.custom_test and run_cpp_comparison_tests
%           are true) 
%     options.cpp_save_to_disk_filename
%         - filename for the hdf5 file saved to disk when using
%           cpp_save_to_disk_only 
%
%     ---------------
%     TEST CASES
%     ---------------
%
%     options.run_source_tests
%         - Boolean to tests with different source types
%     options.run_bin_sensor_tests
%         - Boolean to tests with binary sensor mask
%     options.run_cuboid_corner_tests
%         - Boolean to tests with cuboid corners
%     options.run_cart_sensor_lin_interp_tests
%         - Boolean to tests with a Cartesian sensor mask and linear
%           interpolation (outputs are NOT compared with C++ code)
%     options.run_cart_sensor_nearest_interp_tests
%         - Boolean to tests with a Cartesian sensor mask and nearest
%           neighbour interp (outputs are NOT compared with C++ code)
%     options.run_display_tests
%         - Boolean to tests with different display options (outputs are
%           NOT compared with C++ code) 
%     options.run_time_reversal_tests
%         - Boolean to tests using time reversal
%     options.test_type
%         - integer controlling the test type, where
%               1: Even grid size, PML inside
%               2: Even grid size, PML outside
%               3: Odd grid size, PML inside
%               4: Odd grid size, PML outside
%               5: Custom grid size and PML options, using:
%                   .Nx
%                   .Ny
%                   .Nz
%                   .pml_x_size
%                   .pml_y_size
%                   .pml_z_size
%                   .pml_inside
%     options.test_dim 
%         - integer controlling dimension
%               1: kspaceFirstOrder1D
%               2: kspaceFirstOrder2D
%               3: kspaceFirstOrder3D
%               4: kspaceFirstOrderAS
%     options.data_cast 
%         - data format for MATLAB simulations
%               'off'
%               'single'
%               'gpuArray-single'
%               'gpuArray-double'
%     options.test_case_start_index  
%         - integer to set the start index to skip tests on restart (set to
%           1 to run all tests) 
%
%     ---------------
%     CUSTOM TEST
%     ---------------
%
%     options.custom_test
%         - Boolean controlling whether a custom test should be performed
%     options.custom_test_case
%         - 25 x 1 integer array with the custom test case to run (see test
%           matrix below for documentation)
%
% ABOUT:
%       author      - Bradley Treeby and Jiri Jaros
%       date        - 4th September 2012
%       last update - 10th August 2020
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2020 Bradley Treeby and Jiri Jaros

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

% ASCII generated using http://patorjk.com/software/taag/

%#ok<*UNRCH>
%#ok<*NASGU>

% force diary off
diary off;

% check for k-Wave Toolbox, and add from default location if not present
if exist('getkWavePath', 'file') == 0
    disp('Warning: The k-Wave toolbox is not currently on the MATLAB path');
    disp('         Attempting to add ''/Software/k-Wave'' to the path...');
    try
        addpath('/Software/k-Wave');
    catch ME
        disp('         Not successful!');
        rethrow(ME);
    end
    disp('         Completed successfully!');
end

% get start time
start_time = datetime('now');

% get PC information
computer_info = getComputerInfo;
v = ver('matlab');

% =========================================================================
% SET GRID SIZE
% =========================================================================

% set total number of grid points and the size of the perfectly matched
% layer (PML) 
switch options.test_type
    case 1
        
        % even grid size, pml inside
        NX = 128;           % [grid points]
        NY = 64;            % [grid points]
        NZ = 32;            % [grid points]

        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 10;    % [grid points]
        
        PML_INSIDE = true;  % location of the PML
        
    case 2

        % even grid size, pml outside
        NX = 128;           % [grid points]
        NY = 64;            % [grid points]
        NZ = 32;            % [grid points]

        PML_X_SIZE = 20;    % [grid points]
        if options.test_dim == 4
            PML_Y_SIZE = 16;% [grid points]
        else
            PML_Y_SIZE = 10;% [grid points]
        end
        PML_Z_SIZE = 11;    % [grid points]   
        
        PML_INSIDE = false; % location of the PML
        
    case 3

        % odd grid size, pml inside
        NX = 125;           % [grid points]
        NY = 63;            % [grid points]
        NZ = 35;            % [grid points]
        
        PML_X_SIZE = 20;    % [grid points]
        PML_Y_SIZE = 10;    % [grid points]
        PML_Z_SIZE = 10;    % [grid points]
        
        PML_INSIDE = true;  % location of the PML
        
    case 4

        % odd grid size, pml outside
        NX = 115;           % [grid points]
        NY = 61;            % [grid points]
        NZ = 29;            % [grid points]
        
        PML_X_SIZE = 10;    % [grid points]
        if options.test_dim == 4
            PML_Y_SIZE = 14;% [grid points]
        else
            PML_Y_SIZE = 10;% [grid points]
        end
        PML_Z_SIZE = 10;    % [grid points]
        
        PML_INSIDE = false; % location of the PML
        
    case 5
        
        % custom grid size and pml position
        NX = options.Nx; 	% [grid points]
        NY = options.Ny; 	% [grid points]
        NZ = options.Nz;   	% [grid points]
        
        PML_X_SIZE = options.pml_x_size;    % [grid points]
        PML_Y_SIZE = options.pml_y_size;    % [grid points]
        PML_Z_SIZE = options.pml_z_size;    % [grid points]
        
        PML_INSIDE = options.pml_inside;    % location of the PML

end

% =========================================================================
% SIMULATION LITERALS
% =========================================================================

% define the properties of the propagation medium
C0      = 1540; 	% [m/s]
RHO0    = 1000;     % [kg/m^3]
ALPHA0  = 0.75;     % [dB/(MHz^Y cm)]
if options.test_dim == 4
    Y   = 2;     	% only y = 2 is supported in AS code
else
    Y   = 1.5;     	% power law exponent
end
BONA    = 6;        % parameter of nonlinearity
SC      = 1.3;    	% scale factor for heterogeneity

% define the error tolerance
if options.test_dim == 4
    ERROR_TOLERANCE = 2e-4;
elseif options.use_gpu_code
    ERROR_TOLERANCE = 2e-5;
else
    ERROR_TOLERANCE = 1e-5;
end

% CFL
CFL = 0.1;

% =========================================================================
% CREATE TEST LOG
% =========================================================================

if options.save_test_log 
    date_string = getDateString;
    diary_filename = [options.output_folder 'TESTLOG_' computer_info.computer_name '_' computer_info.operating_system_type '_MATLAB_' v.Version '_' v.Release '_' num2str(options.test_dim) 'D_' options.data_cast '_TESTTYPE_' num2str(options.test_type) '_' date_string];
    diary([diary_filename '.txt']);
end

% =========================================================================
% LIST OF SIMULATION OPTIONS
% =========================================================================

% -----------------------------------------------------------------------
% PARAMETER 01 (LIN)
%   0: Linear
%   1: Nonlinear
% PARAMETER 02 (ABS)
%   0: Lossless
%   1: Power law absorption (y ~= 2)
%   2: Stokes absorption (y = 2)
% PARAMETER 03 (HET)
%   0: Homogeneous
%   1: Heterogeneous c0 and rho0 only
%   2: Heterogeneous
% PARAMETER 04 (SMH)
%   0: 'Smooth' not set (defaults to [true, false, false])
%   1: 'Smooth', true
% PARAMETER 05 (ALP) 
%   0: medium.alpha_mode not set
%   1: medium.alpha_mode = 'no_absorption'
%   2: medium.alpha_mode = 'no_dispersion'
% -----------------------------------------------------------------------
% PARAMETER 06 (SRC)
%   0: initial pressure
%   1: pressure source
%   2: velocity-x source
%   3: velocity-y source (2D and 3D only)
%   4: velocity-z source (3D only)
%   5: velocity-x/y/z source
%   6: transducer source (3D only)
% PARAMETER 07 (MNY)
%   0: single source
%   1: source many
% PARAMETER 08 (DIR)
%   0: additive-no-correction source condition
%   1: dirichlet source condition
%   2: additive source condition
% -----------------------------------------------------------------------
% PARAMETER 09 (BIN)
%   0: binary sensor mask
%   1: Cartesian sensor mask, linear interpolation
%   2: Cartesian sensor mask, nearest neighbour interpolation
%   3: binary sensor mask with sensor directivity (2D only)
%   4: cuboid corners sensor mask
% PARAMETER 10 (P-U)
%   0: record only pressure (no sensor.record input)
%   1: record pressure and velocity
%   2: record max, min and rms pressure and velocity
%   3: record everything BUT u_non_staggered or intensity
%   4: record everything
% PARAMETER 11 (STR) (3D only)
%   0: 'StreamToDisk', false
%   1: 'StreamToDisk', true
%   2: 'StreamToDisk', 50
% PARAMETER 12 (RCT)
%   0: 'DataRecast', false
%   1: 'DataRecast', true
% PARAMETER 13 (FRQ)
%   0: sensor.frequency_response not set
%   1: sensor.frequency_response = [0.5e6, 50]
% PARAMETER 14 (RSI)
%   0: sensor.record_start_index not set
%   1: sensor.record_start_index = 200;
% -----------------------------------------------------------------------
% PARAMETER 15 (LOG)
%   0: 'CreateLog', false (default)
%   1: 'CreateLog', true
% PARAMETER 16 (PML)
%   0: 'PlotPML', true (default)
%   1: 'PlotPML', false (default)
% PARAMETER 17 (LSC)
%   0: 'LogScale', false (default)
%   1: 'LogScale', true
% PARAMETER 18 (DIS)
%   0: 'DisplayMask' not set (defaults to sensor.mask)
%   1: 'DisplayMask', 'off'
%   2: 'DisplayMask', makeBall
% PARAMETER 19 (PFQ)
%   0: 'PlotFreq' not set (defaults to 10) 
%   1: 'PlotFreq', 30
% PARAMETER 20 (PSC)
%   0: 'PlotScale' set to -[p0, p0]
%   1: 'PlotScale' set to 'auto'
% PARAMETER 21 (PLT)
%   0: 'PlotSim' not set (defaults to true) 
%   1: 'PlotSim', false
%   2: 'MeshPlot', true (2D only)
% PARAMETER 22 (LAY)
%   0: 'PlotLayout' not set (defaults to false)
%   1: 'PlotLayout', true
% PARAMETER 23 (MOV)
%   0: 'RecordMovie' not set (defaults to false)
%   1: 'RecordMovie', true
% -----------------------------------------------------------------------
% PARAMETER 24 (TRV)
%   0: run only forward simulation
%   1: run time reversal simulation
%   2: run time reversal including attenuation compensation
% -----------------------------------------------------------------------
% PARAMETER 25 (CPP) (2D and 3D only)
%   0: run MATLAB version only
%   1: run C++ version and compare outputs
% -----------------------------------------------------------------------

% set the number of variations for each parameters
variations = [2, 3, 3, 2, 3,...
    7, 2, 2,...
    5, 5, 3, 2, 2,...
    2, 2, 2, 3, 2, 2, 2, 2, 3,...
    3,...
    2];

% =========================================================================
% SET TEST CASES
% =========================================================================

% set pseudonym for running C++ comparisons - for tests that make sense to
% compare with the C++ code, option 25 is set to C, otherwise option 25 is
% set to 0
options.run_cpp_comparison_tests = (options.run_cpp_comparison_tests) && (options.test_dim ~= 1);
C = double(options.run_cpp_comparison_tests);

test_cases = [];
if options.run_source_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   A: TEST SOURCE CONDITIONS IN HOMOGENEOUS AND HETEROGENEOUS MEDIA
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   6   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...      
    0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    1   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   2   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   2   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...         
    1   0   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   2   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   2   2   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   0   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    0   1   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   1   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...   
    0   2   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...     
    1   2   0   0   0   3   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...       
% -----------------------------------------------------------------------------------------------------
    ];
end

if options.run_corrected_source_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   A: TEST SOURCE CONDITIONS IN HOMOGENEOUS AND HETEROGENEOUS MEDIA
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   1   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   1   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   2   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   3   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   4   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   5   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   1   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   2   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   3   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   4   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   1   0   0   5   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
% -----------------------------------------------------------------------------------------------------
    ];
end

if options.run_bin_sensor_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   B: TEST BINARY SENSOR CONDITIONS
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...      
    0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   0   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...    
    0   0   0   0   0   0   0   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   0   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   0   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   0   4   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
% -----------------------------------------------------------------------------------------------------
    ];
end        
    
if options.run_cart_sensor_lin_interp_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   C: TEST CARTESIAN SENSOR CONDITIONS WITH LINEAR INTERPOLATION
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   1   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   1   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   1   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   1   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   1   4   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   1   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...   
% -----------------------------------------------------------------------------------------------------
    ];
end 

if options.run_cart_sensor_nearest_interp_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   D: TEST CARTESIAN SENSOR CONDITIONS WITH NEAREST NEIGHBOUR INTERPOLATION
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   2   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   2   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   2   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   2   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   2   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   2   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   2   4   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   2   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0; ...  
% -----------------------------------------------------------------------------------------------------
    ];
end 

if options.run_cuboid_corner_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   E: TEST CUBOID CORNERS SENSOR MASK
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   0   0   0   4   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   4   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   4   3   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   4   4   0   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   4   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...      
    0   0   0   0   0   0   0   0   4   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   4   2   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...    
    0   0   0   0   0   0   0   0   4   3   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   4   4   0   1   0   0   0   0   0   0   0   0   0   0   0   0   C; ...    
    0   0   0   0   0   0   0   0   4   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   4   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   4   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0; ...       
    0   0   0   0   0   0   0   0   4   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0; ...           
    0   0   0   0   0   0   0   0   4   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   4   1   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
    0   0   0   0   0   0   0   0   4   2   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ... 
    0   0   0   0   0   0   0   0   4   3   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   4   4   0   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...
    0   0   0   0   0   0   0   0   4   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   C; ...  
% -----------------------------------------------------------------------------------------------------
    ];
end 

if options.run_display_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   F: TEST DISPLAY CONDITIONS
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0   0   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0   1   0   0   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   1   0   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   1   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   1   0   0; ...
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0   0   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0   0   0   0; ...      
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0; ...    
    0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0; ...  
    0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0; ...  
    0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0; ...  
    0   2   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0; ...      
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0; ...  
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0; ...      
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0; ...    
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0   1   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0   1   0   0; ... 
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0; ...     
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0; ...         
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   2   0   0   0   0; ...         
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   2   0   0   0   0; ...             
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   2   0   0   0   0; ...             
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   1   0   0; ...         
% -----------------------------------------------------------------------------------------------------
    ]; 
end 

if options.run_time_reversal_tests
    test_cases = [test_cases;...
% -----------------------------------------------------------------------------------------------------
%   G: TEST TIME REVERSAL
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% -----------------------------------------------------------------------------------------------------
    0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   C; ...     
    0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ... 
    0   0   0   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ...   
    0   0   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   C; ...     
    0   0   1   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ... 
    0   0   1   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   1   0; ...       
    0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   1   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   1   0   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...   
    0   1   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   1   1   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   1   1   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...
    0   2   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   2   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   2   0   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...   
    0   2   1   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...     
    0   2   1   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ... 
    0   2   1   0   0   0   0   0   2   0   0   1   0   0   0   0   0   0   0   0   0   0   0   2   0; ...       
% -----------------------------------------------------------------------------------------------------
    ]; 
end 

% =========================================================================
% REMOVE TESTS THAT DON'T MAKE SENSE FOR PARTICULAR DIMENSIONS
% =========================================================================

% remove power law absorption tests for axisymmetric
if options.test_dim == 4
    test_cases(find(test_cases(:, 2) == 1), :) = []; %#ok<*FNDSB>
end

% remove transducer source tests
if options.test_dim ~= 3
    test_cases(find(test_cases(:, 6) == 6), :) = [];
end

% remove velocity-x/yz source tests
if options.test_dim == 1
    test_cases(find(test_cases(:, 6) == 5), :) = [];
end

% remove velocity-z source tests
if options.test_dim ~= 3
    test_cases(find(test_cases(:, 6) == 4), :) = [];
end

% remove velocity-y source tests
if options.test_dim == 1
    test_cases(find(test_cases(:, 6) == 3), :) = [];
end

% remove directivity tests
if options.test_dim ~= 2 || strcmp(options.data_cast, 'gpuArray-single') || strcmp(options.data_cast, 'gpuArray-double') 
    test_cases(find(test_cases(:, 9) == 3), :) = [];
end

% remove Cart interp nearest tests without recast
if ~(strcmp(options.data_cast, 'off') || strcmp(options.data_cast, 'single'))
    test_cases(find( (test_cases(:, 9) == 2) & (test_cases(:, 12) == 0)), :) = [];
end

% remove Cart sensor mask test (interp1 not supported by parallel computing
% toolbox)
if (options.test_dim == 1) && (strcmp(options.data_cast, 'gpuArray-double') || strcmp(options.data_cast, 'gpuArray-single'))
    test_cases(find(test_cases(:, 9) == 1), :) = [];
end

% remove stream to disk tests
if options.test_dim ~= 3
    test_cases(find(test_cases(:, 11) ~= 0), :) = [];
end

% remove frequency response tests for GPU tests without recasting
if ~(strcmp(options.data_cast, 'off') || strcmp(options.data_cast, 'single'))
    test_cases(find( (test_cases(:, 13) == 1) & (test_cases(:, 12) == 0)), :) = [];
end

% remove tests that record intensity with the GPU without recasting
if ~(strcmp(options.data_cast, 'off') || strcmp(options.data_cast, 'single'))
    test_cases(find( (test_cases(:, 10) == 4) & (test_cases(:, 12) == 0)), :) = [];
end

% remove mesh plot tests
if ~((options.test_dim == 2) || (options.test_dim == 4))
    test_cases(find(test_cases(:, 21) == 2), :) = [];
end

% remove plot layout tests
if strcmp(options.data_cast, 'GPUsingle') || strcmp(options.data_cast, 'GPUdouble')
    test_cases(find(test_cases(:, 22) == 1), :) = [];
end

% =========================================================================
% REPLACE WITH A SINGLE TEST CASE IF DESIRED
% =========================================================================

if options.custom_test
    test_cases = options.custom_test_case; 
end

% =========================================================================
% COMMAND LINE UPDATE
% =========================================================================
clc;
disp('   ');
disp('-----------------------------------------------------------------------------------------------');
disp('                 _      __        __                _____         _            ');
disp('                | | __  \ \      / /_ ___   _____  |_   _|__  ___| |_ ___ _ __ ');
disp('                | |/ /___\ \ /\ / / _` \ \ / / _ \   | |/ _ \/ __| __/ _ \ ''__|');
disp('                |   <_____\ V  V / (_| |\ V /  __/   | |  __/\__ \ ||  __/ |   ');
disp('                |_|\_\     \_/\_/ \__,_| \_/ \___|   |_|\___||___/\__\___|_|   ');
disp('  ');                                                                
disp('-----------------------------------------------------------------------------------------------');
disp('  ');
disp(['START DATE:                                      ' char(start_time)]);
disp(['COMPUTER:                                        ' computer_info.computer_name]);
disp(['USER NAME:                                       ' computer_info.user_name]);
disp(['O/S:                                             ' computer_info.operating_system]);
disp(['O/S TYPE:                                        ' computer_info.operating_system_type]);
disp(['MATLAB VERSION:                                  ' v.Version ' ' v.Release]);
disp(['K-WAVE VERSION:                                  ' computer_info.kwave_version]);
disp('  ');
disp(['DATA CAST:                                       ' options.data_cast]);
disp(['TEST DIM:                                        ' num2str(options.test_dim)]);
switch options.test_dim
    case 1
disp(['GRID SIZE:                                       ' num2str(NX)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE)]);
    case 2
disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE)]);
    case 3
disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY) ' by ' num2str(NZ)]);
disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE) ' by ' num2str(PML_Z_SIZE)]);
end
disp(['PML INSIDE:                                      ' num2str(PML_INSIDE)]);
disp(['NONUNIFORM GRID:                                 ' num2str(options.use_nonuniform_grid)]);
disp('  ');
disp(['RUN SOURCE TESTS:                                ' num2str(options.run_source_tests)]);
disp(['RUN BININARY SENSOR TESTS:                       ' num2str(options.run_bin_sensor_tests)]);
disp(['RUN CARTESIAN SENSOR TESTS (LIN INTERP):         ' num2str(options.run_cart_sensor_lin_interp_tests)]);
disp(['RUN CARTESIAN SENSOR TESTS (NN INTERP):          ' num2str(options.run_cart_sensor_nearest_interp_tests)]);
disp(['RUN DISPLAY TESTS:                               ' num2str(options.run_display_tests)]);
disp(['RUN TIME REVERSAL TESTS:                         ' num2str(options.run_time_reversal_tests)]);
disp(['RUN C++ COMPARISON TESTS:                        ' num2str(options.run_cpp_comparison_tests)]);
disp('  ');
disp('  ');
disp(['The total number of possible test combinations is: ' num2str(prod(variations))]);
disp(['This would take an estimated ' scaleTime(prod(variations)*30) ' to test at 30 seconds per simulation']);
disp(['The number of tested combinations is: ' num2str(size(test_cases, 1))]);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% reset the total number of errors
number_errors = 0;
number_cpp_errors = 0;
location_errors = [];
location_cpp_errors = [];

% extract parameters from test cast settings
for comp_index = options.test_case_start_index:size(test_cases, 1)

    disp('  ');
    disp('---------------------------------------------------------------------------------------------------');
    disp(['INDEX OF CURRENT TEST = ' num2str(comp_index)])
    disp('---------------------------------------------------------------------------------------------------');
    disp('LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP');
    disp('01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25');
    disp('---------------------------------------------------------------------------------------------------');
    disp(num2str(test_cases(comp_index, :), '%1d   '));
    disp('---------------------------------------------------------------------------------------------------');
    disp(' ');
    
    % -------------------------------------------------------------------
    
    % PARAMETER 1
    NONLINEAR = (test_cases(comp_index, 1) == 1);
    
    % PARAMETER 2
    ABSORBING = test_cases(comp_index, 2);

    % PARAMETER 3
    if test_cases(comp_index, 3) == 0
        HETEROGENEOUS_RHO0   = false;
        HETEROGENEOUS_C0     = false;
        HETEROGENEOUS_BONA   = false;    % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = false;    % not used if ABSORBING = false
    elseif test_cases(comp_index, 3) == 1
        HETEROGENEOUS_RHO0   = true;
        HETEROGENEOUS_C0     = true;
        HETEROGENEOUS_BONA   = false;    % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = false;    % not used if ABSORBING = false
    elseif test_cases(comp_index, 3) == 2
        HETEROGENEOUS_RHO0  = true;
        HETEROGENEOUS_C0     = true;
        HETEROGENEOUS_BONA   = true;     % not used if NONLINEAR = false
        HETEROGENEOUS_ALPHA0 = true;     % not used if ABSORBING = false
    else
        error('unknown setting');
    end

    % PARAMETER 4
    SMOOTH = (test_cases(comp_index, 4) == 1);
    
    % PARAMETER 5
    switch test_cases(comp_index, 5)
        case 0
            ALPHA_MODE = [];
        case 1
            ALPHA_MODE = 'no_absorption';
        case 2
            ALPHA_MODE = 'no_dispersion';
    end    
    
	% -------------------------------------------------------------------
    
    % PARAMETER 6
    SOURCE_TYPE = test_cases(comp_index, 6);

    % PARAMETER 7
    SOURCE_MANY = test_cases(comp_index, 7);
    
    % PARAMETER 8
    SOURCE_MODE = test_cases(comp_index, 8);
    
    % -------------------------------------------------------------------
        
    % PARAMETER 9
    switch test_cases(comp_index, 9)
        case 0
            BINARY_SENSOR_MASK = true;
            CUBOID_CORNERS     = false;            
            SENSOR_DIRECTIVITY = false;
        case 1
            BINARY_SENSOR_MASK = false;      
            CUBOID_CORNERS     = false;             
            CART_INTERP = 'linear';
        case 2
            BINARY_SENSOR_MASK = false;    
            CUBOID_CORNERS     = false;             
            CART_INTERP = 'nearest';
        case 3
            BINARY_SENSOR_MASK = true;
            CUBOID_CORNERS     = false;            
            SENSOR_DIRECTIVITY = true;
        case 4
            BINARY_SENSOR_MASK = false;
            CUBOID_CORNERS     = true;
    end
    
    % PARAMETER 10
    switch test_cases(comp_index, 10)
        case 0
            SENSOR_RECORD_FIELDS = [];
            T_END = 10e-6;
        case 1
            SENSOR_RECORD_FIELDS = {'p', 'u'};
            T_END = 10e-6;
        case 2
            SENSOR_RECORD_FIELDS = {'p_max', 'p_rms', 'p_min', 'u_max', 'u_rms', 'u_min'};
            T_END = 10e-6;
        case 3
            SENSOR_RECORD_FIELDS = {'p', 'p_max', 'p_rms', 'p_min', 'p_max_all', 'p_min_all', 'p_final', ...
                                    'u', 'u_max', 'u_rms', 'u_min', 'u_max_all', 'u_min_all', 'u_final'};
            T_END = 10e-6;
        case 4
            SENSOR_RECORD_FIELDS = {'p', 'p_max', 'p_rms', 'p_min', 'p_max_all', 'p_min_all', 'p_final', ...
                                    'u', 'u_max', 'u_rms', 'u_min', 'u_max_all', 'u_min_all', 'u_final',...
                                    'u_non_staggered', 'I', 'I_avg'};
            T_END = 10e-6;
    end
    
    % PARAMETER 11
    switch test_cases(comp_index, 11)
        case 0
            STREAM_TO_DISK = [];
        case 1
            STREAM_TO_DISK = true;
        case 2
            STREAM_TO_DISK = 20;
    end
    
    % PARAMETER 12
    DATA_RECAST = (test_cases(comp_index, 12) == 1);
        
    % PARAMETER 13
    switch test_cases(comp_index, 13)
        case 0
            FREQUENCY_RESPONSE = [];
        case 1
            FREQUENCY_RESPONSE = [0.5e6, 50];
    end
    
    % PARAMETER 14
    switch test_cases(comp_index, 14)
        case 0
            RECORD_START_INDEX = 1;
        case 1
            RECORD_START_INDEX = 200;
    end
    
    % -------------------------------------------------------------------
    
    % PARAMETER 15
    switch test_cases(comp_index, 15)
        case 0
            CREATE_LOG = [];
        case 1
            CREATE_LOG = true;
    end
    
    % PARAMETER 16
    switch test_cases(comp_index, 16)
        case 0
            PLOT_PML = [];
        case 1
            PLOT_PML = false;
    end
    
    % PARAMETER 17
    switch test_cases(comp_index, 17)
        case 0
            PLOT_LOG_SCALE = [];
        case 1
            PLOT_LOG_SCALE = false;
    end
    
    % PARAMETER 18
    switch test_cases(comp_index, 18)
        case 0
            DISPLAY_MASK = [];
        case 1
            DISPLAY_MASK = 'off';
        case 2
            DISPLAY_MASK = 'ball';
    end
    
    % PARAMETER 19
    switch test_cases(comp_index, 19)
        case 0
            PLOT_FREQ = [];
        case 1
            PLOT_FREQ = 30;
    end    
    
    % PARAMETER 20
    switch test_cases(comp_index, 20)
        case 0
            PLOT_SCALE = [];
        case 1
            PLOT_SCALE = 'auto';
    end 
    
    % PARAMETER 21
    switch test_cases(comp_index, 21)
        case 0
            PLOT_SIM = [];
            MESH_PLOT = [];
        case 1
            PLOT_SIM = false;
            MESH_PLOT = [];
        case 2
            PLOT_SIM = [];
            MESH_PLOT = true;
            
    end 
    
    % PARAMETER 22
    switch test_cases(comp_index, 22)
        case 0
            PLOT_LAYOUT = [];
        case 1
            PLOT_LAYOUT = true;
    end     
        
    % PARAMETER 23
    switch test_cases(comp_index, 23)
        case 0
            SAVE_MOVIE = [];
        case 1
            SAVE_MOVIE = true;
    end  
    
    % -------------------------------------------------------------------
    
    % PARAMETER 24
    RUN_TIME_REVERSAL = test_cases(comp_index, 24);
    
    % -------------------------------------------------------------------
    
    % PARAMETER 25
    if options.test_dim > 1
        RUN_CPP_CODE = test_cases(comp_index, 25);
    else
        RUN_CPP_CODE = 0;
    end
       
    % -------------------------------------------------------------------
    
    % run simulation
    try
        run_simulation;
    catch ME
        if options.continue_past_errors
            disp(' ');
            disp('  SIMULATION ERROR:');
            disp('  -----------------');
            disp(['  message identifier: ' ME.identifier]);
            disp(['  occured in: ' ME.stack(1).name ' at line ' num2str(ME.stack(1).line)]);
            disp(['  message: ' ME.message]);
            disp(' ');
            number_errors = number_errors + 1;
            location_errors = [location_errors, comp_index]; %#ok<AGROW>
        else
            diary off;
            rethrow(ME);
        end
    end
    
    % split test log in parts (some problems with really long log files)
    if options.save_test_log && ~rem(comp_index, options.split_log_after_n_tests)
        diary off;
        diary([diary_filename '_PART' num2str(floor(comp_index / options.split_log_after_n_tests) + 1) '.txt']);
    end

end

disp(' ');
disp('---------------------------------------------------------------------------------------------------');
disp(['TOTAL RUNTIME ERRORS:        ' num2str(number_errors)])
disp(['INDEX OF RUNTIME ERRORS:     ' num2str(location_errors)]);
if options.run_cpp_comparison_tests
disp(['C++ SIMS WITH L_INF > ' num2str(ERROR_TOLERANCE) ':  ' num2str(number_cpp_errors)]);
disp(['INDEX OF C++ ERRORS:         ' num2str(location_cpp_errors)]);
end
disp(['ELAPSED TIME:                ' scaleTime(seconds(datetime('now') - start_time))]);
disp('---------------------------------------------------------------------------------------------------');

% switch off log
if options.save_test_log
    
    % close diary
    diary off;
    pause(1);
    diary off;
    pause(1);
    diary off;
    
end

% save summary
if options.save_test_log
    
    % open new diary file
    diary([diary_filename '_summary.txt']);
    
    % print summary
    disp('-----------------------------------------------------------------------------------------------');
    disp('                 _      __        __                _____         _            ');
    disp('                | | __  \ \      / /_ ___   _____  |_   _|__  ___| |_ ___ _ __ ');
    disp('                | |/ /___\ \ /\ / / _` \ \ / / _ \   | |/ _ \/ __| __/ _ \ ''__|');
    disp('                |   <_____\ V  V / (_| |\ V /  __/   | |  __/\__ \ ||  __/ |   ');
    disp('                |_|\_\     \_/\_/ \__,_| \_/ \___|   |_|\___||___/\__\___|_|   ');
    disp('  ');                                                                
    disp('-----------------------------------------------------------------------------------------------');
    disp('  ');
    disp(['START DATE:                                      ' char(start_time)]);
    disp(['COMPUTER:                                        ' computer_info.computer_name]);
    disp(['USER NAME:                                       ' computer_info.user_name]);
    disp(['O/S:                                             ' computer_info.operating_system]);
    disp(['O/S TYPE:                                        ' computer_info.operating_system_type]);
    disp(['MATLAB VERSION:                                  ' v.Version ' ' v.Release]);
    disp(['K-WAVE VERSION:                                  ' computer_info.kwave_version]);
    disp('  ');
    disp(['DATA CAST:                                       ' options.data_cast]);
    disp(['TEST DIM:                                        ' num2str(options.test_dim)]);
    switch options.test_dim
        case 1
    disp(['GRID SIZE:                                       ' num2str(NX)]);
    disp(['PML SIZE:                                        ' num2str(PML_X_SIZE)]);
        case 2
    disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY)]);
    disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE)]);
        case 3
    disp(['GRID SIZE:                                       ' num2str(NX) ' by ' num2str(NY) ' by ' num2str(NZ)]);
    disp(['PML SIZE:                                        ' num2str(PML_X_SIZE) ' by ' num2str(PML_Y_SIZE) ' by ' num2str(PML_Z_SIZE)]);
    end
    disp(['PML INSIDE:                                      ' num2str(PML_INSIDE)]);
    disp(['NONUNIFORM GRID:                                 ' num2str(options.use_nonuniform_grid)]);
    disp('  ');
    disp(['RUN SOURCE TESTS:                                ' num2str(options.run_source_tests)]);
    disp(['RUN BININARY SENSOR TESTS:                       ' num2str(options.run_bin_sensor_tests)]);
    disp(['RUN CARTESIAN SENSOR TESTS (LIN INTERP):         ' num2str(options.run_cart_sensor_lin_interp_tests)]);
    disp(['RUN CARTESIAN SENSOR TESTS (NN INTERP):          ' num2str(options.run_cart_sensor_nearest_interp_tests)]);
    disp(['RUN DISPLAY TESTS:                               ' num2str(options.run_display_tests)]);
    disp(['RUN TIME REVERSAL TESTS:                         ' num2str(options.run_time_reversal_tests)]);
    disp(['RUN C++ COMPARISON TESTS:                        ' num2str(options.run_cpp_comparison_tests)]);
    disp('  ');
    disp(['NUMBER OF TESTS:                                 ' num2str(size(test_cases, 1))]);
    disp('  ');
    disp('---------------------------------------------------------------------------------------------------');
    disp(['TOTAL RUNTIME ERRORS:        ' num2str(number_errors)])
    disp(['INDEX OF RUNTIME ERRORS:     ' num2str(location_errors)]);
    if options.run_cpp_comparison_tests
    disp(['C++ SIMS WITH L_INF > ' num2str(ERROR_TOLERANCE) ':  ' num2str(number_cpp_errors)]);
    disp(['INDEX OF C++ ERRORS:         ' num2str(location_cpp_errors)]);
    end
    disp(['ELAPSED TIME:                ' scaleTime(seconds(datetime('now') - start_time))]);
    disp('---------------------------------------------------------------------------------------------------');
    
    % close diary file
    diary off;
    
end

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% ========================================================================
% ------------------------------------------------------------------------
% ========================================================================
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

function run_simulation()
% Nested subfunction to set options and run simulations. For 3D
% simulations, the outputs from the MATLAB and C++ version of the code are
% compared. (Note, nested subfunctions can see all variables created
% above.) 

% create empty input arguments;
input_args = {};

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% create heterogeneous region
switch options.test_dim
    case 1
        heterog = ones(NX, 1);
        heterog(round(NX/2) - 9:round(NX/2) + 10) = 1;
    case 2
        heterog = makeDisc(NX, NY, round(NX/2), round(NY/2), 10);
    case 3
        heterog = makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), 10);
    case 4
        heterog = makeDisc(NX, 2 * NY, round(NX/2), round(NY/2), 10);
        heterog = heterog(:, NY + 1:end);
end

% set density
if HETEROGENEOUS_RHO0
    medium.density = RHO0*ones(size(heterog));
    medium.density(heterog == 1) = RHO0*SC;
else
    medium.density = RHO0; 
end

% set sound speed
if HETEROGENEOUS_C0
    medium.sound_speed = C0*ones(size(heterog));
    medium.sound_speed(heterog == 1) = C0*SC;
else
    medium.sound_speed = C0;
end
    
% set BONA
if NONLINEAR
    if HETEROGENEOUS_BONA
        medium.BonA = BONA*ones(size(heterog));
        medium.BonA(heterog == 1) = BONA*SC;
    else
        medium.BonA = BONA;
    end
end

% set absorption terms
if ABSORBING
    
    % assign medium properties
    if HETEROGENEOUS_ALPHA0
        medium.alpha_coeff  = ALPHA0*ones(size(heterog));
        medium.alpha_coeff(heterog == 1) = ALPHA0*SC;
        medium.alpha_power = Y;
    else
        medium.alpha_coeff = ALPHA0;      
        medium.alpha_power = Y;
    end
    
    % set absorption mode
    if ~isempty(ALPHA_MODE)
        medium.alpha_mode = ALPHA_MODE;
    end
    
    % set stokes absorption
    if ABSORBING == 2
        medium.alpha_mode = 'stokes';
        medium.alpha_power = 2;
    end
    
end

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/NX;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
switch options.test_dim
    case 1
        kgrid = kWaveGrid(NX, dx);
    case {2,4}
        kgrid = kWaveGrid(NX, dx, NY, dy);
    case 3
        kgrid = kWaveGrid(NX, dx, NY, dy, NZ, dz);
end

% =========================================================================
% DEFINE NONUNIFORM GRID SETTINGS
% =========================================================================

% if options.use_nonuniform_grid
%     
%    ...
%     
% end

% =========================================================================
% DEFINE THE TIME ARRAY
% =========================================================================

% create the time array
kgrid.makeTime(medium.sound_speed, CFL, T_END);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa] (this is also defined as a literal in the errorNorms function!)
tone_burst_freq = 0.5e6;    	% [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
input_signal = source_strength * toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% =========================================================================
% DEFINE THE SOURCE CONDITIONS
% =========================================================================

% 0: initial pressure
% 1: pressure source
% 2: velocity-x source
% 3: velocity-y source
% 4: velocity-z source
% 5: velocity-x/y/z source
% 6: transducer source

switch SOURCE_TYPE
    case 0
        
        % -----------------------------
        % OPTION 0: Initial Pressure
        % ----------------------------- 
        
        % create ball shaped initial pressure distribution
        source_radius = 5;
        switch options.test_dim
            case 1
                source.p0 = ones(NX, 1);
                source.p0(round(NX/2) - source_radius + 1:round(NX/2) + source_radius) = source_strength;
            case {2,4}
                source.p0 = source_strength * makeDisc(NX, NY, round(NX/2), round(NY/2), source_radius);
            case 3
                source.p0 = source_strength * makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), source_radius);
        end
        
    case 1
        
        % -----------------------------
        % OPTION 1: Pressure Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch options.test_dim
            case 1
                source.p_mask = zeros(NX, 1);
                source.p_mask(PML_X_SIZE + 1) = 1;
                source.p_mask(PML_X_SIZE + source_radius) = 1;
            case 2
                source.p_mask = zeros(NX, NY);
                source.p_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius) = 1;
            case 3
                source.p_mask = zeros(NX, NY, NZ);
                source.p_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, round(NZ/2) - round(source_radius/2) + 1:round(NZ/2) + round(source_radius/2)) = 1;
            case 4
                source.p_mask = zeros(NX, NY);
                source.p_mask(PML_X_SIZE + 1, 1:source_radius * 3) = 1;
        end

        % assign the source term
        if SOURCE_MANY
            switch options.test_dim
                case 1
                    focus_position = 0;
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
                case {2,4}
                    focus_position = [-5*dx, 0];
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.p = focus(kgrid, input_signal, source.p_mask, focus_position, C0);
            end
        else
            source.p = input_signal;
        end
        
    case 2
        
        % -----------------------------
        % OPTION 2: Velocity X Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch options.test_dim
            case 1
                source.u_mask = zeros(NX, 1);
                source.u_mask(PML_X_SIZE + 1) = 1;
                source.u_mask(PML_X_SIZE + source_radius) = 1;
            case 2
                source.u_mask = zeros(NX, NY);
                source.u_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius) = 1;
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(PML_X_SIZE + 1, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, round(NZ/2) - source_radius/2 + 1:round(NZ/2) + source_radius/2) = 1;
            case 4
                source.u_mask = zeros(NX, NY);
                source.u_mask(PML_X_SIZE + 1, 1:source_radius * 3) = 1;                
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch options.test_dim
                case 1
                    focus_position = 0;
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case {2,4}
                    focus_position = [-5*dx, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 3
                    focus_position = [0, 0, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
            end
        else
            source.ux = input_signal./(C0*RHO0);
        end
        
    case 3
        
        % -----------------------------
        % OPTION 3: Velocity Y Source
        % -----------------------------
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch options.test_dim
            case 2
                source.u_mask = zeros(NX, NY);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, PML_Y_SIZE + 1) = 1;
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, PML_Y_SIZE + 1, round(NZ/2) - round(source_radius/2) + 1:round(NZ/2) + round(source_radius/2)) = 1;
            case 4
                source.u_mask = zeros(NX, NY);
                source.u_mask(PML_X_SIZE + 1:PML_X_SIZE + source_radius * 2, PML_Y_SIZE + 1) = 1;
            otherwise
                error('source.uy not supported in 1D');
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch options.test_dim
                case 2
                    focus_position = [-5*dx, 0];
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 3
                    focus_position = [-5*dx, 0, 0];
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                case 4
                    focus_position = [5*dx, 5*dx];
                    source.uy = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
            end
        else
            source.uy = input_signal./(C0*RHO0);
        end
        
    case 4
        
        % -----------------------------
        % OPTION 4: Velocity Z Source
        % -----------------------------        
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch options.test_dim
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, PML_Z_SIZE + 1) = 1;
            otherwise
                error('source.uz only supported in 3D');
        end

        % assign the source term scaled by the impedance
        if SOURCE_MANY
            focus_position = [0, 0, 0];
            source.uz = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
        else
            source.uz = input_signal./(C0*RHO0);
        end
        
    case 5
        
        % -----------------------------
        % OPTION 5: Velocity XYZ Source
        % ----------------------------- 
        
        % set the source mask to be a small rectangle
        source_radius = 15;
        switch options.test_dim
            case 2
                source.u_mask = zeros(NX, NY);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, PML_Y_SIZE + 1) = 1;
            case 3
                source.u_mask = zeros(NX, NY, NZ);
                source.u_mask(round(NX/2) - source_radius + 1:round(NX/2) + source_radius, round(NY/2) - source_radius + 1:round(NY/2) + source_radius, PML_Z_SIZE + 1) = 1;
            case 4
                source.u_mask = zeros(NX, NY);
                source.u_mask(PML_X_SIZE + 1:PML_X_SIZE + source_radius * 2, PML_Y_SIZE + 1) = 1;
            otherwise
                error('x/y/z velocity source not supported in 1D');
        end
        
        % assign the source term scaled by the impedance
        if SOURCE_MANY
            switch options.test_dim
                case {2,4}
                    focus_position = [-5*dx, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                    source.uy = source.ux;
                case 3
                    focus_position = [0, 0, 0];
                    source.ux = focus(kgrid, input_signal./(C0*RHO0), source.u_mask, focus_position, C0);
                    source.uy = source.ux;
                    source.uz = source.ux;
            end
        else
            switch options.test_dim
                case {2,4}
                    source.ux = input_signal./(C0*RHO0);
                    source.uy = input_signal./(C0*RHO0);
                case 3
                    source.ux = input_signal./(C0*RHO0);
                    source.uy = input_signal./(C0*RHO0);
                    source.uz = input_signal./(C0*RHO0);
            end
        end
        
    case 6
        
        % -----------------------------
        % OPTION 6: Transducer Source
        % -----------------------------
        
        % scale the source magnitude by the source_strength divided by the
        % impedance (the source is assigned to the particle velocity)
        input_signal = input_signal./(C0*RHO0);

        % physical properties of the transducer
        transducer.number_elements = 32;    % total number of transducer elements
        transducer.element_width = 1;       % width of each element [grid points]
        transducer.element_length = 12;     % length of each element [grid points]
        transducer.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points]
        transducer.radius = inf;            % radius of curvature of the transducer [m]

        % calculate the width of the transducer in grid points
        transducer_width = transducer.number_elements*transducer.element_width ...
            + (transducer.number_elements - 1)*transducer.element_spacing;

        % use this to position the transducer in the middle of the computational grid
        transducer.position = round([PML_X_SIZE + 1, round(NY/2) - transducer_width/2, round(NZ/2) - transducer.element_length/2]);

        % properties used to derive the beamforming delays
        transducer.sound_speed = 1540;              % sound speed [m/s]
        transducer.focus_distance = 20e-3;          % focus distance [m]
        transducer.elevation_focus_distance = 19e-3;% focus distance in the elevation plane [m]
        transducer.steering_angle = 0;              % steering angle [degrees]

        % apodization
        transducer.transmit_apodization = 'Rectangular';    
        transducer.receive_apodization = 'Rectangular';

        % define the transducer elements that are currently active
        transducer.active_elements = ones(transducer.number_elements, 1);

        % append input signal used to drive the transducer
        transducer.input_signal = input_signal;

        % create the transducer using the defined settings
        transducer = kWaveTransducer(kgrid, transducer);
        
        % pass the transducer to the source
        source = transducer;
        
end

% source condition for pressure sources
if (SOURCE_TYPE == 1)
    switch SOURCE_MODE
        case 0
            source.p_mode = 'additive-no-correction';
        case 1
            source.p_mode = 'dirichlet';
        case 2
            source.p_mode = 'additive';
    end
end
    
% source condition for pressure sources
if (SOURCE_TYPE > 1) && (SOURCE_TYPE < 6)
    switch SOURCE_MODE
        case 0
            source.u_mode = 'additive-no-correction';
        case 1
            source.u_mode = 'dirichlet';
        case 2
            source.u_mode = 'additive';
    end
end

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

if BINARY_SENSOR_MASK

    % define a sensor mask through the central plane, offset by a small
    % amount
    switch options.test_dim
        case 1
            sensor.mask = zeros(NX, 1);
            sensor.mask(round(NX/2):round(NX/2)+1) = 1;
        case 2
            sensor.mask = zeros(NX, NY);
            sensor.mask(:, round(NY/2) + round(NY/8)) = 1;
            
            if SENSOR_DIRECTIVITY
                sensor.directivity_angle = zeros(NX, NY);
                sensor.directivity_angle(sensor.mask == 1) = 0;
                sensor.directivity_size = 3*kgrid.dx;
            end
        case 3
            sensor.mask = zeros(NX, NY, NZ);
            sensor.mask(:, :, round(NZ/2) + round(NZ/8)) = 1;
        case 4
            sensor.mask = zeros(NX, NY);
            sensor.mask(:, round(NY/2) + round(NY/8), :) = 1;
    end

elseif CUBOID_CORNERS
    
    % define two intersecting cuboid corners
    switch options.test_dim
        case 1
            sensor.mask = [round(NX/4), round(NX/2 + NX/8);...
                           round(NX/2), round(3*NX/4)].';
        case {2,4}
            sensor.mask = [round(NX/4), round(NY/4), round(NX/2 + NX/8), round(NY/2);...
                           round(NX/2), round(NY/4 + NY/8), round(3*NX/4), round(3*NY/4)].';
        case 3
            sensor.mask = [round(NX/4), round(NY/4), round(NZ/4), round(NX/2 + NX/8), round(NY/2), round(NZ/2);...
                           round(NX/2), round(NY/4 + NY/8), round(NZ/4), round(3*NX/4), round(3*NY/4), round(NZ/2)].';
    end    
    
else
    
    % define a series of Cartesian points to collect the data
    switch options.test_dim
        case 1
            x = (-22:2:22)*dx;          % [m]
            sensor.mask = x;
        case 2
            x = (-22:2:22)*dx;          % [m]
            y = (-11:1:11)*dy;          % [m]
            sensor.mask = [x; y];
        case 3
            x = (-22:2:22)*dx;          % [m]
            y = (-11:1:11)*dy;          % [m]
            z = 10*dz*ones(size(x));    % [m]
            sensor.mask = [x; y; z];
        case 4
            x = (-22:2:22)*dx;          % [m]
            y = (0:1:22)*dy;            % [m]
            sensor.mask = [x; y];
    end
    
    % add interpolation option to input arguments
    input_args = [input_args, {'CartInterp', CART_INTERP}];
    
end

% sensor record option
if ~isempty(SENSOR_RECORD_FIELDS)
    sensor.record = SENSOR_RECORD_FIELDS;
end

% sensor frequency response
if ~isempty(FREQUENCY_RESPONSE)
    sensor.frequency_response = FREQUENCY_RESPONSE;
end

% sensor record start time
if RECORD_START_INDEX > 1
    sensor.record_start_index = RECORD_START_INDEX;
end

% =========================================================================
% SET OPTIONAL INPUT ARGUMENTS
% =========================================================================

% stream to disk option
if ~isempty(STREAM_TO_DISK)
    input_args = [input_args, {'StreamToDisk', STREAM_TO_DISK}];
end

% save logfile
if ~isempty(CREATE_LOG)
    input_args = [input_args, {'CreateLog', CREATE_LOG}];
end

% plot pml
if ~isempty(PLOT_PML)
    input_args = [input_args, {'PlotPML', PLOT_PML}];
end

% logscale plot
if ~isempty(PLOT_LOG_SCALE)
    input_args = [input_args, {'LogScale', PLOT_LOG_SCALE}];
end

% display mask
if ~isempty(DISPLAY_MASK)
    if strcmp(DISPLAY_MASK, 'ball')
        switch options.test_dim
            case 1
                DISPLAY_MASK = zeros(NX, 1);
                DISPLAY_MASK(round(NX/2):round(NX/2)+1) = 1;
            case {2,4}
                DISPLAY_MASK = makeDisc(NX, NY, round(NX/2), round(NY/2), 10);
            case 3
                DISPLAY_MASK = makeBall(NX, NY, NZ, round(NX/2), round(NY/2), round(NZ/2), 10);
        end
    end
    input_args = [input_args, {'DisplayMask', DISPLAY_MASK}];
end

% plot frequency
if ~isempty(PLOT_FREQ)
    input_args = [input_args, {'PlotFreq', PLOT_FREQ}];
end

% plot scale
if isempty(PLOT_SCALE)
    if options.test_dim == 1
        PLOT_SCALE = [-source_strength*1.2, source_strength*1.2];
    else
        PLOT_SCALE = [-source_strength/2, source_strength/2];
    end
end
    
% plot on or off
if options.force_plot_off
    input_args = [input_args, {'PlotSim', false}];
elseif ~isempty(PLOT_SIM)
    input_args = [input_args, {'PlotSim', PLOT_SIM}];
end
    
% plot layout
if ~isempty(PLOT_LAYOUT) && ~options.force_plot_off 
    input_args = [input_args, {'PlotLayout', PLOT_LAYOUT}];
end

% save movie
if ~isempty(SAVE_MOVIE) && ~options.force_plot_off 
    movie_name = [options.output_folder getDateString '_kWave_Test_' num2str(comp_index)];
    input_args = [input_args, {'RecordMovie', SAVE_MOVIE, 'MovieName', movie_name}];
end

% smoothing
if SMOOTH
    input_args = [input_args, {'Smooth', SMOOTH}];
end

% mesh plot
if ~isempty(MESH_PLOT)
    input_args = [input_args, {'MeshPlot', true}];
end

% PML size
switch options.test_dim
    case 1
        input_args = [input_args, {'PMLSize', PML_X_SIZE}];
    case {2,4}
        input_args = [input_args, {'PMLSize', [PML_X_SIZE, PML_Y_SIZE]}];
    case 3
        input_args = [input_args, {'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE]}];
end

% set additional input settings
input_args = [input_args, {'PMLInside', PML_INSIDE,...
    'DataCast', options.data_cast, 'DataRecast', DATA_RECAST, 'PlotScale', PLOT_SCALE}];

% set additional input settings for C++ codes
input_args_cpp = {...
    'NumThreads', options.number_threads, ...
    'VerboseLevel', options.cpp_verbose_level, ...
    'ThreadBinding', options.cpp_thread_binding, ...
    };
if options.use_gpu_code
    input_args_cpp = [input_args_cpp, {'BinaryName', options.cuda_binary_name, 'BinaryPath', options.cuda_binary_path}];
else
    input_args_cpp = [input_args_cpp, {'BinaryName', options.cpp_binary_name, 'BinaryPath', options.cpp_binary_path}];
end

% =========================================================================
% RUN THE MATLAB SIMULATION
% =========================================================================

% skip MATLAB simulation if options.custom_test and
% options.cpp_save_to_disk_only are true

if ~(options.custom_test && options.cpp_save_to_disk_only)
    
    % run the simulation using k-Wave
    switch options.test_dim
        case 1
            sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
        case 2
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        case 3
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
        case 4
            sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:});
    end

    % run the time reversal simulation 
    if RUN_TIME_REVERSAL

        % assign the time reversal data
        sensor_time_rev = sensor;
        if isstruct(sensor_data)
            sensor_time_rev.time_reversal_boundary_data = sensor_data.p;
        else
            sensor_time_rev.time_reversal_boundary_data = sensor_data;
        end

        % attenuation compensation options
        if RUN_TIME_REVERSAL == 2

            % define the cutoff frequency for the filter
            f_cutoff = 1e6;

            % create the filter to regularise the absorption parameters
            if ~(options.test_dim == 4 || ABSORBING == 2)
                medium.alpha_filter = getAlphaFilter(kgrid, medium, f_cutoff);
            end

            % reverse the sign of the absorption proportionality coefficient
            medium.alpha_sign = [-1, 1];        % [absorption, dispersion];

        end

        % run the reconstruction
        switch options.test_dim
            case 1
                p0_recon = kspaceFirstOrder1D(kgrid, medium, [], sensor_time_rev, input_args{:});
            case 2
                p0_recon = kspaceFirstOrder2D(kgrid, medium, [], sensor_time_rev, input_args{:});
            case 3
                p0_recon = kspaceFirstOrder3D(kgrid, medium, [], sensor_time_rev, input_args{:});
            case 4
                p0_recon = kspaceFirstOrderAS(kgrid, medium, [], sensor_time_rev, input_args{:});
        end

        % clean up
        clear sensor_time_rev

    end
end

% =========================================================================
% RUN THE C++ SIMULATION
% =========================================================================

% compare the output of the MATLAB code with the C++ if required
if RUN_CPP_CODE 
    
    % save the input file to disk and stop if required
    if (options.custom_test && options.cpp_save_to_disk_only)
        switch options.test_dim
            case 2
                kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', options.cpp_save_to_disk_filename);
            case 3
                kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', options.cpp_save_to_disk_filename);
            case 4
                kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', options.cpp_save_to_disk_filename);
            otherwise
                error('Save to disk option not supported in 1D.');
        end
        return
    end
    
    % run the simulation again using C++ version
    switch options.test_dim
        case 2
            if options.use_gpu_code
                sensor_data_cpp = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            else
                sensor_data_cpp = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            end
        case 3
            if options.use_gpu_code
                sensor_data_cpp = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            else
                sensor_data_cpp = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            end
        case 4
            if options.use_gpu_code
                sensor_data_cpp = kspaceFirstOrderASG(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            else
                sensor_data_cpp = kspaceFirstOrderASC(kgrid, medium, source, sensor, input_args{:}, input_args_cpp{:});
            end
        otherwise
            error('Save to disk option not supported in 1D.');
    end
    
    if RUN_TIME_REVERSAL

        % assign the time reversal data from the C++ forward simulation
        sensor_time_rev = sensor;
        if isstruct(sensor_data_cpp)
            sensor_time_rev.time_reversal_boundary_data = sensor_data_cpp.p;
        else
            sensor_time_rev.time_reversal_boundary_data = sensor_data_cpp;
        end

        % run the reconstruction
        switch options.test_dim
            case 2
                if options.use_gpu_code
                    p0_recon_cpp = kspaceFirstOrder2DG(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                else
                    p0_recon_cpp = kspaceFirstOrder2DC(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                end                
            case 3
                if options.use_gpu_code
                    p0_recon_cpp = kspaceFirstOrder3DG(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                else
                    p0_recon_cpp = kspaceFirstOrder3DC(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                end
            case 4
                if options.use_gpu_code
                    p0_recon_cpp = kspaceFirstOrderASG(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                else
                    p0_recon_cpp = kspaceFirstOrderASC(kgrid, medium, [], sensor_time_rev, input_args{:}, input_args_cpp{:});
                end
        end
    end

    % set plot axis
    x_axis = (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3;
    y_axis = kgrid.y_vec*1e3;
    
    % number of time steps
    if isfield(sensor, 'record_start_index')
        NT = kgrid.Nt - sensor.record_start_index + 1;
    else
        NT = kgrid.Nt;
    end
    
    % check to see if data should be recast to the CPU
    if ismember(options.data_cast, {'gsingle', 'GPUsingle', 'gpuArray-single'})
        recast = @(x) single(x);
    else
        recast = @(x) x;
    end
    
    % set command line output
    disp(' ');
    disp(' ');
    disp('  C++ ACCURACY COMPARED TO MATLAB:');
    disp('  --------------------------------');
    
    % compare outputs for p0_recon
    if RUN_TIME_REVERSAL
        
        % assign the output data (it shouldn't need reshaping)
        mat = p0_recon;
        cpp = p0_recon_cpp;
        
        % display error norms
        [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'p0_recon', comp_index, number_cpp_errors, location_cpp_errors);
        
        % plot
        if options.plot_cpp_errors
            mat = squeeze(mat(round(end/2), :, :));
            cpp = squeeze(cpp(round(end/2), :, :));
            plot_title = 'p0 recon';
            plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
        end
        
    else

        % set the number of cuboids (set to 1 if cuboids aren't used)
        if CUBOID_CORNERS
            cuboid_num = 2;
        else
            cuboid_num = 1;
        end
        
        % loop over output data for each cuboid
        for cuboid_index = 1:cuboid_num
        
            % compare the outputs for p
            if (isfield(sensor, 'record') && ismember('p', sensor.record)) || ~isfield(sensor, 'record')

                % reshape the output data
                if ~CUBOID_CORNERS
                    if isfield(sensor, 'record')
                        if (kgrid.dim == 2)
                            mat = reshape(recast(sensor_data(cuboid_index).p), [NX, NT]);
                            cpp = reshape(recast(sensor_data_cpp(cuboid_index).p), [NX, NT]);                          
                        else                              
                            mat = reshape(recast(sensor_data(cuboid_index).p), [NX, NY, NT]);                          
                            cpp = reshape(recast(sensor_data_cpp(cuboid_index).p), [NX, NY, NT]);
                        end
                    else
                        if (kgrid.dim == 2)
                            mat = reshape(recast(sensor_data), [NX, NT]);
                            cpp = reshape(recast(sensor_data_cpp), [NX, NT]);
                        else                            
                            mat = reshape(recast(sensor_data), [NX, NY, NT]);
                            cpp = reshape(recast(sensor_data_cpp), [NX, NY, NT]);
                        end
                    end
                else                    
                    mat = recast(sensor_data(cuboid_index).p);
                    cpp = recast(sensor_data_cpp(cuboid_index).p);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').p'], comp_index, number_cpp_errors, location_cpp_errors);

                % calculate the distribution of fundamental and harmonic
                if NONLINEAR && options.plot_cpp_errors && ~CUBOID_CORNERS

                    % compute the frequency axis
                    freq = (0:ceil((NT + 1)/2) - 1) ./ (kgrid.dt*NT);

                    % compute the index at which the source frequency and its harmonics occur
                    [~, f0_index] = findClosest(freq, tone_burst_freq);
                    [~, f1_index] = findClosest(freq, tone_burst_freq*2);

                    % preallocate the beam pattern variables
                    beam_pattern_f0 = zeros(NX, NY);
                    beam_pattern_f1 = zeros(NX, NY);
                    beam_pattern_total = zeros(NX, NY);
                    beam_pattern_f0_cpp = zeros(NX, NY);
                    beam_pattern_f1_cpp = zeros(NX, NY);
                    beam_pattern_total_cpp = zeros(NX, NY);

                    % compute the amplitude spectrum of the time series recorded at each sensor
                    % point, and then extract the corresponding amplitudes at the fundamental
                    % frequency and second harmonic.
                    for x_index = 1:NX
                        for y_index = 1:NY

                            % compute the amplitude spectrum
                            amp_spect = spect(squeeze(mat(x_index, y_index, :)), 1/kgrid.dt); %, 'Window', 'Hanning');
                            amp_spect_cpp = spect(squeeze(cpp(x_index, y_index, :)), 1/kgrid.dt);

                            % extract the amplitude at the source frequency and store
                            beam_pattern_f0(x_index, y_index) = amp_spect(f0_index);
                            beam_pattern_f0_cpp(x_index, y_index) = amp_spect_cpp(f0_index);

                            % extract the amplitude at the source frequency and store
                            beam_pattern_f1(x_index, y_index) = amp_spect(f1_index); 
                            beam_pattern_f1_cpp(x_index, y_index) = amp_spect_cpp(f1_index);

                            % extract the integral of the total amplitude spectrum
                            beam_pattern_total(x_index, y_index) = sum(amp_spect(:));
                            beam_pattern_total_cpp(x_index, y_index) = sum(amp_spect_cpp(:));

                        end
                    end

                    % produce plots
                    mat = beam_pattern_total;
                    cpp = beam_pattern_total_cpp;
                    plot_title = 'Total Pressure';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);

                    mat = beam_pattern_f0;
                    cpp = beam_pattern_f0_cpp;          
                    plot_title = 'Fundamental';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);

                    mat = beam_pattern_f1;
                    cpp = beam_pattern_f1_cpp; 
                    plot_title = '2nd Harmonic';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);

                elseif options.plot_cpp_errors

                    % compute the maximum value to plot
                    plot_title = 'Total Pressure';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);

                end
            end

            % compare the outputs for p_max
            if isfield(sensor, 'record') && ismember('p_max', sensor.record)

                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).p_max);
                        cpp = recast(sensor_data_cpp(cuboid_index).p_max);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).p_max), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).p_max), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).p_max);
                    cpp = recast(sensor_data_cpp(cuboid_index).p_max);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').p_max'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'p max';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

            end

            % compare the outputs for p_min
            if isfield(sensor, 'record') && ismember('p_min', sensor.record)

                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2                        
                        mat = recast(sensor_data(cuboid_index).p_min);
                        cpp = recast(sensor_data_cpp(cuboid_index).p_min);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).p_min), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).p_min), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).p_min);
                    cpp = recast(sensor_data_cpp(cuboid_index).p_min);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').p_min'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'p min';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

            end        

            % compare the outputs for p_rms
            if isfield(sensor, 'record') && ismember('p_rms', sensor.record)

                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2                        
                        mat = recast(sensor_data(cuboid_index).p_rms);
                        cpp = recast(sensor_data_cpp(cuboid_index).p_rms);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).p_rms), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).p_rms), [NX, NY]);
                    end
                    
                else
                    mat = recast(sensor_data(cuboid_index).p_rms);
                    cpp = recast(sensor_data_cpp(cuboid_index).p_rms);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').p_rms'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'p rms';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

            end    

            % compare outputs for u
            if isfield(sensor, 'record') && ismember('u', sensor.record)

                % ------
                % UX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = reshape(recast(sensor_data(cuboid_index).ux), [NX,  NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux), [NX, NT]);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).ux), [NX, NY, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux), [NX, NY, NT]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).ux);
                    cpp = recast(sensor_data_cpp(cuboid_index).ux);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').ux'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'ux';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = reshape(recast(sensor_data(cuboid_index).uy), [NX, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy), [NX, NT]);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).uy), [NX, NY, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy), [NX, NY, NT]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).uy);
                    cpp = recast(sensor_data_cpp(cuboid_index).uy);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uy'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uy';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UZ
                % ------
                
                if kgrid.dim == 3
                
                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).uz), [NX, NY, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uz), [NX, NY, NT]);
                    else
                        mat = recast(sensor_data(cuboid_index).uz);
                        cpp = recast(sensor_data_cpp(cuboid_index).uz);
                    end

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uz'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'uz';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end
                    
                end

            end

            % compare outputs for u non staggered
            if isfield(sensor, 'record') && ismember('u_non_staggered', sensor.record)

                % ------
                % UX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).ux_non_staggered);
                        cpp = recast(sensor_data_cpp(cuboid_index).ux_non_staggered);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).ux_non_staggered), [NX, NY, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux_non_staggered), [NX, NY, NT]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).ux_non_staggered);
                    cpp = recast(sensor_data_cpp(cuboid_index).ux_non_staggered);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').ux_non_staggered'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'ux_non_staggered';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                   if kgrid.dim == 2
                     mat = recast(sensor_data(cuboid_index).uy_non_staggered);
                     cpp = recast(sensor_data_cpp(cuboid_index).uy_non_staggered);
                   else
                      mat = reshape(recast(sensor_data(cuboid_index).uy_non_staggered), [NX, NY, NT]);
                      cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy_non_staggered), [NX, NY, NT]);
                   end
                else
                    mat = recast(sensor_data(cuboid_index).uy_non_staggered);
                    cpp = recast(sensor_data_cpp(cuboid_index).uy_non_staggered);
                end

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uy_non_staggered'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uy_non_staggered';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UZ
                % ------
                
                if kgrid.dim == 3

                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).uz_non_staggered), [NX, NY, NT]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uz_non_staggered), [NX, NY, NT]);
                    else
                          mat = recast(sensor_data(cuboid_index).uz_non_staggered);
                          cpp = recast(sensor_data_cpp(cuboid_index).uz_non_staggered);
                    end

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uz_non_staggered'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'uz_non_staggered';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end
                    
                end

            end
            
            % compare the outputs for u_max
            if isfield(sensor, 'record') && ismember('u_max', sensor.record)

                % ------
                % UX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).ux_max);
                        cpp = recast(sensor_data_cpp(cuboid_index).ux_max);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).ux_max), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux_max), [NX, NY]);
                    end
                        
                else
                    mat = recast(sensor_data(cuboid_index).ux_max);
                    cpp = recast(sensor_data_cpp(cuboid_index).ux_max);
                end                    

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').ux_max'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'ux max';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).uy_max);
                        cpp = recast(sensor_data_cpp(cuboid_index).uy_max);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).uy_max), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy_max), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).uy_max);
                    cpp = recast(sensor_data_cpp(cuboid_index).uy_max);
                end   

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uy_max'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uy max';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UZ
                % ------
                
                if kgrid.dim == 3    
                
                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).uz_max), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uz_max), [NX, NY]);
                    else
                        mat = recast(sensor_data(cuboid_index).uz_max);
                        cpp = recast(sensor_data_cpp(cuboid_index).uz_max);
                    end   

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uz_max'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'uz max';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end       
                    
                end

            end

            % compare the outputs for u_min
            if isfield(sensor, 'record') && ismember('u_min', sensor.record)

                % ------
                % UX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).ux_min);
                        cpp = recast(sensor_data_cpp(cuboid_index).ux_min);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).ux_min), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux_min), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).ux_min);
                    cpp = recast(sensor_data_cpp(cuboid_index).ux_min);
                end                       

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').ux_min'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'ux min';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).uy_min);
                        cpp = recast(sensor_data_cpp(cuboid_index).uy_min);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).uy_min), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy_min), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).uy_min);
                    cpp = recast(sensor_data_cpp(cuboid_index).uy_min);
                end 

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uy_min'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uy min';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                
                % ------
                % UZ
                % ------
                
                if kgrid.dim == 3
                
                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).uz_min), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uz_min), [NX, NY]);
                    else
                        mat = recast(sensor_data(cuboid_index).uz_min);
                        cpp = recast(sensor_data_cpp(cuboid_index).uz_min);
                    end 

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uz_min'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'uz min';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end
                end

            end        

            % compare the outputs for u_rms
            if isfield(sensor, 'record') && ismember('u_rms', sensor.record)

                % ------
                % UX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).ux_rms);
                        cpp = recast(sensor_data_cpp(cuboid_index).ux_rms);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).ux_rms), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).ux_rms), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).ux_rms);
                    cpp = recast(sensor_data_cpp(cuboid_index).ux_rms);
                end 

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').ux_rms'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'ux rms';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).uy_rms);
                        cpp = recast(sensor_data_cpp(cuboid_index).uy_rms);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).uy_rms), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uy_rms), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).uy_rms);
                    cpp = recast(sensor_data_cpp(cuboid_index).uy_rms);
                end 

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uy_rms'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uy rms';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % UZ
                % ------
                
                if kgrid.dim == 3
                
                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).uz_rms), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).uz_rms), [NX, NY]);
                    else
                        mat = recast(sensor_data(cuboid_index).uz_rms);
                        cpp = recast(sensor_data_cpp(cuboid_index).uz_rms);
                    end 

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').uz_rms'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'uz rms';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end
                    
                end

            end    

            % compare outputs for I_avg
            if isfield(sensor, 'record') && ismember('I_avg', sensor.record)

                % ------
                % IX
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).Ix_avg);
                        cpp = recast(sensor_data_cpp(cuboid_index).Ix_avg);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).Ix_avg), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).Ix_avg), [NX, NY]);
                    end
                else
                    mat = recast(sensor_data(cuboid_index).Ix_avg);
                    cpp = recast(sensor_data_cpp(cuboid_index).Ix_avg);
                end 

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').Ix_avg'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'Ix avg';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % IY
                % ------
                
                % reshape the output data
                if ~CUBOID_CORNERS
                    if kgrid.dim == 2
                        mat = recast(sensor_data(cuboid_index).Iy_avg);
                        cpp = recast(sensor_data_cpp(cuboid_index).Iy_avg);
                    else
                        mat = reshape(recast(sensor_data(cuboid_index).Iy_avg), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).Iy_avg), [NX, NY]);
                    end 
                else
                    mat = recast(sensor_data(cuboid_index).Iy_avg);
                    cpp = recast(sensor_data_cpp(cuboid_index).Iy_avg);
                end 

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').Iy_avg'], comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'Iy avg';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end

                % ------
                % IZ
                % ------
                
                if kgrid.dim == 3
                
                    % reshape the output data
                    if ~CUBOID_CORNERS
                        mat = reshape(recast(sensor_data(cuboid_index).Iz_avg), [NX, NY]);
                        cpp = reshape(recast(sensor_data_cpp(cuboid_index).Iz_avg), [NX, NY]);
                    else
                        mat = recast(sensor_data(cuboid_index).Iz_avg);
                        cpp = recast(sensor_data_cpp(cuboid_index).Iz_avg);
                    end 

                    % display error norms
                    [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, ['sensor_data(' num2str(cuboid_index) ').Iz_avg'], comp_index, number_cpp_errors, location_cpp_errors);

                    % plot
                    if options.plot_cpp_errors
                        plot_title = 'Iz avg';
                        plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                    end
                    
                end

            end   
        end 
        
        % compare variables that don't depend on cuboids, i.e., capture the
        % whole matrix

        % compare the outputs for p_final
        if isfield(sensor, 'record') && ismember('p_final', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).p_final);
            cpp = recast(sensor_data_cpp(1).p_final);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).p_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'p final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

        end      

        % compare the outputs for p_max_all
        if isfield(sensor, 'record') && ismember('p_max_all', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).p_max_all);
            cpp = recast(sensor_data_cpp(1).p_max_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).p_max_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'p_max_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

        end    

        % compare the outputs for p_min_all
        if isfield(sensor, 'record') && ismember('p_min_all', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).p_min_all);
            cpp = recast(sensor_data_cpp(1).p_min_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).p_min_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'p_min_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

        end           
        
        % compare the outputs for u_final
        if isfield(sensor, 'record') && ismember('u_final', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).ux_final);
            cpp = recast(sensor_data_cpp(1).ux_final);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).ux_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'ux final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

            % reshape the output data
            mat = recast(sensor_data(1).uy_final);
            cpp = recast(sensor_data_cpp(1).uy_final);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uy_final', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'uy final';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

            if kgrid == 3
                
                % reshape the output data   
                mat = recast(sensor_data(1).uz_final);
                cpp = recast(sensor_data_cpp(1).uz_final);

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uz_final', comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uz final';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end
                
            end

        end

        % compare the outputs for u_max_all
        if isfield(sensor, 'record') && ismember('u_max_all', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).ux_max_all);
            cpp = recast(sensor_data_cpp(1).ux_max_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).ux_max_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'ux max_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

            % reshape the output data
            mat = recast(sensor_data(1).uy_max_all);
            cpp = recast(sensor_data_cpp(1).uy_max_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uy_max_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'uy max_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

            if kgrid == 3
                
                % reshape the output data
                mat = recast(sensor_data(1).uz_max_all);
                cpp = recast(sensor_data_cpp(1).uz_max_all);

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uz_max_all', comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uz max_all';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end    
                 
            end

        end        

        % compare the outputs for u_min_all
        if isfield(sensor, 'record') && ismember('u_min_all', sensor.record)

            % reshape the output data
            mat = recast(sensor_data(1).ux_min_all);
            cpp = recast(sensor_data_cpp(1).ux_min_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).ux_min_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'ux min_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end

            % reshape the output data
            mat = recast(sensor_data(1).uy_min_all);
            cpp = recast(sensor_data_cpp(1).uy_min_all);

            % display error norms
            [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uy_min_all', comp_index, number_cpp_errors, location_cpp_errors);

            % plot
            if options.plot_cpp_errors
                plot_title = 'uy min_all';
                plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
            end
            
            if kgrid == 3
                
                % reshape the output data
                mat = recast(sensor_data(1).uz_min_all);
                cpp = recast(sensor_data_cpp(1).uz_min_all);

                % display error norms
                [number_cpp_errors, location_cpp_errors] = errorNorms(mat, cpp, 'sensor_data(1).uz_min_all', comp_index, number_cpp_errors, location_cpp_errors);

                % plot
                if options.plot_cpp_errors
                    plot_title = 'uz min_all';
                    plotErrors(x_axis, y_axis, mat, cpp, comp_index, plot_title, options.save_cpp_comparison_plots_to_disk, [options.output_folder options.image_foldername]);
                end
                
            end

        end         
        
    end
end

% end for subfunction
end

% subfunction to calculate and display error norms
function [num_errors, loc_errors] = errorNorms(mat, cpp, error_title, test_index, num_errors, loc_errors)
   
% check for nans
mat_nans = any(isnan(mat(:)));
cpp_nans = any(isnan(cpp(:)));
   
% get the max values for reference
max_mat = max(abs(mat(:)));
max_cpp = max(abs(cpp(:)));

% compute the errors
L2 = sqrt( sum(abs(mat(:).^2 - cpp(:).^2)) / sum(mat(:).^2) );
LINF = max(abs(mat(:) - cpp(:)));
LINF_REL = LINF/max_mat;

% display the error norms
disp(' ');
disp(['  Error in ' error_title]);
disp(['    MAX VALS = ' num2str(max_mat) ' (MATLAB) ' num2str(max_cpp) ' (CPP)']);
disp(['    L2 = ' num2str(L2)]);
disp(['    LINF = ' num2str(LINF) ' (' num2str(LINF_REL) ' normalised to max value)' ]);

% display nan warnings
if mat_nans
    disp('    NANS detected in MATLAB output');
end
if cpp_nans
    disp('    NANS detected in C++ output');
end

% add error to list
if (min(LINF, LINF_REL) > ERROR_TOLERANCE) || cpp_nans || mat_nans
    num_errors = num_errors + 1;
    if isempty(loc_errors)
        loc_errors = test_index;
    elseif loc_errors(end) ~= test_index
        loc_errors = [loc_errors, test_index];
    end
end

end

% subfunction to produce plots
function plotErrors(x_axis, y_axis, matlab_data, cpp_data, test_index, plot_title, save_plots, image_folder)
    
    % collapse inputs if multi-dimensional
    if numDim(matlab_data) == 4
        
        % reshape [X, Y, Z, T] into [pos_index, time_index]
        matlab_data = reshape(matlab_data, [], size(matlab_data, 4));
        cpp_data    = reshape(cpp_data, [], size(cpp_data, 4));
        
    elseif numDim(matlab_data) == 3
        if size(matlab_data, 1) == NX
            
            % plot central slice
            matlab_data = squeeze(matlab_data(round(end/2), :, :));
            cpp_data    = squeeze(cpp_data(round(end/2), :, :));
            
        else
            
            % take MIP over time
            matlab_data = max(matlab_data, [], 3);
            cpp_data    = max(cpp_data, [], 3);
            
        end
    end

    % produce plot of errors
    h = figure;
    subplot(1, 4, 1);
    imagesc(y_axis, x_axis, matlab_data/1e6);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title(['k-Wave (' plot_title ')']);
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;

    subplot(1, 4, 2);
    imagesc(y_axis, x_axis, cpp_data/1e6);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title(['C++ (' plot_title ')']);
    colormap(jet(256));
    c = colorbar;
    ylabel(c, 'Pressure [MPa]');
    axis image;

    subplot(1, 4, 3);
    imagesc(y_axis, x_axis, 100*abs(matlab_data - cpp_data)./matlab_data);
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title('Local Error');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, '[%]');
    axis image;

    subplot(1, 4, 4);
    imagesc(y_axis, x_axis, 100*abs(matlab_data - cpp_data)./max(abs(matlab_data(:))));
    set(gca, 'FontSize', 10);
    xlabel('y-position [mm]');
    ylabel('x-position [mm]');
    title('Global Error');
    colormap(jet(256));
    c = colorbar;
    ylabel(c, '[%]');
    axis image;
    
    % scale figure
    scaleFig(2, 1);
    
    % save plots to disk and close
    if save_plots
        
        % check folder exists
        if ~exist(image_folder, 'dir')
            mkdir(image_folder);
        end
        
        % save plot
        print(h, '-dpng','-r300', [image_folder 'test' num2str(test_index) '-' plot_title  '.png']);
        close all hidden; 
        
    end
    
end

% end for parent function
end