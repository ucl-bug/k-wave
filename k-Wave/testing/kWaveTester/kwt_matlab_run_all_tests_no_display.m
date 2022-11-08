% Script to call kWavetester to test the functionality of
% kspaceFirstOrder1D, kspaceFirstOrder2D, and kspaceFirstOrder3D with the
% display turned off.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 3rd May 2017

clearvars;

% load defaults
kwt_load_defaults;

% no comparison with C++ code
options.run_cpp_comparison_tests = false;

% no display
options.run_display_tests = false;
options.force_plot_off = true;

% double precision
options.data_cast = 'off'; 
kwt_loop_dim;

% single precision
options.data_cast = 'single'; 
kwt_loop_dim;