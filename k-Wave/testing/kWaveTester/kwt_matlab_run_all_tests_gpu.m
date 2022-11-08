% Script to call kWavetester to test the functionality of
% kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kspaceFirstOrderAS using the 'DataCast' option set to run on the GPU.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 15th April 2019

clearvars;

% load defaults
kwt_load_defaults;

% stop when an error is encountered?
options.continue_past_errors = true;

% no comparison with C++ code
options.run_cpp_comparison_tests = false;

% single precision GPU using parallel computing toolbox
options.data_cast = 'gpuArray-single'; 
kwt_loop_dim;

% double precision GPU using parallel computing toolbox
options.data_cast = 'gpuArray-double'; 
kwt_loop_dim;