% Script to call kWaveTester to test the functionality of
% kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kspaceFirstOrderAS.
%
% author: Bradley Treeby
% date: 29th July 2019
% last update: 29th July 2019

clearvars;

% load defaults
kwt_load_defaults;

% stop when an error is encountered
options.continue_past_errors                    = false;

% turn on comparisons with C++ code
options.run_cpp_comparison_tests                = true;

% turn off display
options.force_plot_off                          = true;

% don't plot differences
options.plot_cpp_errors                         = false;

% turn off tests that don't compare the C++ code
options.run_display_tests                       = false;
options.run_corrected_source_tests              = false;
options.run_cart_sensor_lin_interp_tests        = false;
options.run_cart_sensor_nearest_interp_tests    = false;

% data cast option for test comparison
options.data_cast                               = 'single';  

% 2D CPU
options.use_gpu_code                            = false;
options.test_dim                                = 2;
kwt_loop_test_type;

% 2D GPU
options.use_gpu_code                            = true;
options.test_dim                                = 2;
kwt_loop_test_type;

% 3D CPU
options.use_gpu_code                            = false;
options.test_dim                                = 3;
kwt_loop_test_type;

% 3D GPU
options.use_gpu_code                            = true;
options.test_dim                                = 3;
kwt_loop_test_type;

% Axisymmetric CPU
options.use_gpu_code                            = false;
options.test_dim                                = 4;
kwt_loop_test_type;