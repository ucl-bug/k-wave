% Script to call kWavetester to test the functionality of
% kspaceFirstOrder-CUDA compared with kspaceFirstOrder2D.
%
% author: Bradley Treeby
% date: 6th March 2019
% last update: 6th March 2019

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

% use the gpu code
options.use_gpu_code                            = true;

% set folder for images (set plot_cpp_errors = true)
options.image_foldername                        = ['cuda_comparison_tests_' getDateString filesep];

% data cast option for test comparison
options.data_cast                               = 'single';  

% set grid dimension
options.test_dim                                = 2;

% run tests
kwt_loop_test_type;