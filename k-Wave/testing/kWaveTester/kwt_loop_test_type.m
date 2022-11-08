% Script to call kWaveTester for all test types.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 23rd February 2017

% 1: Even grid size, PML inside
options.test_type = 1;
kWaveTester(options);
close all hidden;

% 2: Even grid size, PML outside
options.test_type = 2;
kWaveTester(options);
close all hidden;

% 3: Odd grid size, PML inside
options.test_type = 3;
kWaveTester(options);
close all hidden;

% 4: Odd grid size, PML outside
options.test_type = 4;
kWaveTester(options);
close all hidden;