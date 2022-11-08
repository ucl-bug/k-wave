% Script to use kWaveTester to run a single custom test case.
%
% author: Bradley Treeby
% date: 3rd March 2019
% last update: 3rd March 2019

clearvars;

% load defaults
kwt_load_defaults;

% turn off log
options.save_test_log = false;

% stop for errors
options.continue_past_errors = false;

% options for plotting
options.force_plot_off = false;
options.plot_cpp_errors = true;
options.save_cpp_comparison_plots_to_disk = false;

% set C++ options
options.number_threads = 'all';
options.cpp_verbose_level = 2;
options.cpp_thread_binding = 1;

% set dimension
options.test_dim = 4;

% test type
options.test_type = 1;

% set this to true to run a single test as defined below
options.custom_test = true;

% custom test options (only used if options.custom_test = true)
% -----------------------------------------------------------------------------------------------------
options.custom_test_case = [...
% -----------------------------------------------------------------------------------------------------
%   LIN ABS HET SMH ALP SRC MNY DIR BIN P-U STR RCT FRQ RSI LOG PML LSC DIS PFQ PSC PLT LAY MOV TRV CPP   
%   01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
% ----------------------------------------------------------------------------------------------------- 
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1; ...  
% -----------------------------------------------------------------------------------------------------
];

% -----------------------------------------------------------------------
% PARAMETER 01 (LIN)
%   0: Linear
%   1: Nonlinear
% PARAMETER 02 (ABS)
%   0: Lossless
%   1: Absorbing
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
%   0: additive source condition
%   1: dirichlet source condition
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
%   2: 'RecordMovie', true / 'MovieType', 'image' (2D only)
% -----------------------------------------------------------------------
% PARAMETER 24 (TRV)
%   0: run only forward simulation
%   1: run time reversal simulation
%   2: run time reversal including attenuation compensation
% -----------------------------------------------------------------------
% PARAMETER 25 (CPP) (3D only)
%   0: run MATLAB version only
%   1: run C++ version and compare outputs
% -----------------------------------------------------------------------

% run test
kWaveTester(options);