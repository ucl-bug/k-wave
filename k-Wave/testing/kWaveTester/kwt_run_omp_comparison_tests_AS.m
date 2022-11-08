% Script to call kWavetester to test the functionality of
% kspaceFirstOrder-OMP compared with kspaceFirstOrderAS.
%
% author: Bradley Treeby
% date: 15th October 2018
% last update: 10th July 2019

clearvars;

% load defaults
kwt_load_defaults;

% no display
options.run_display_tests                       = false;
options.force_plot_off                          = true;

% option to plot the differences
options.plot_cpp_errors                         = false;

% turn off tests that don't compare the C++ code
options.run_corrected_source_tests              = false;
options.run_cart_sensor_lin_interp_tests        = false;
options.run_cart_sensor_nearest_interp_tests    = false;

% stop if errors are encountered
options.continue_past_errors                    = false;

% set folder for images (set plot_cpp_errors = true)
options.image_foldername                        = ['omp_comparison_tests_AS_' getDateString filesep];

% set grid dimension
options.test_dim                                = 4;

% datacast option for MATLAB code
options.data_cast                               = 'single';

% run tests
kwt_loop_test_type;