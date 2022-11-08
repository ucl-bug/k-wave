% Script to call kwt_loop_test_type for all dimensions
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 16th April 2018

% 1D: kspaceFirstOrder1D
options.test_dim = 1;
kwt_loop_test_type;

% 2D: kspaceFirstOrder2D
options.test_dim = 2;
kwt_loop_test_type;

% 3D: kspaceFirstOrder3D
options.test_dim = 3;
kwt_loop_test_type;

% AS: kspaceFirstOrderAS
options.test_dim = 4;
kwt_loop_test_type;