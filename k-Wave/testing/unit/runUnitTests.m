function test_struct = runUnitTests(wildcard)
%RUNUNITTESTS Run MATLAB unit tests.
%
% DESCRIPTION:
%     runUnitTests runs all unit tests (i.e., all .m files) within the
%     folder k-Wave/testing/unit/. The unit tests are expected to be
%     in the form:
%
%         pass = unitTest(plot_simulations, plot_comparisons)
%
%     where
%  
%         pass             - Boolean denoting whether the test passed
%         plot_simulations - Boolean controlling whether any simulations
%                            are plotted
%         plot_comparisons - Boolean controlling whether any comparisons
%                            (e.g., error norms) are plotted
%
%     After all unit tests have been performed, a summary of the test names
%     and pass status is displayed.
%
% OPTIONAL INPUTS:
%         wildcard         - String with wildcard pattern to match test
%                            filenames
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 17th February 2014
%     last update   - 4th March 2025
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2025 Bradley Treeby

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

% start the clock
regression_start_time = clock;

% literals
plot_simulations = 'false';
plot_comparisons = 'false';

% get a list of m-files in the directory
filenames = what;
filenames = filenames.m;

% remove any files that start with 'runUnitTests'
filenames(startsWith(filenames, 'runUnitTests')) = [];

% filter filenames based on wildcard
if nargin > 0
    filenames = filenames(contains(filenames, wildcard));
end

filenames=filenames(1:2)
% extract number of files to test
num_files = length(filenames);

% keep a list of whether the test passed or failed
test_result = false(num_files, 1);

% =========================================================================
% RUN TESTS
% =========================================================================

% run the examples one by one
for filename_index = 1:num_files

    % remove test pass variable
    clear test_pass;
    
    % trim the .m extension
    fn = filenames{filename_index};
    fn = fn(1:end - 2);
    
    % display the filename
    disp(['Running ' fn ' (Test ' num2str(filename_index) ' of ' num2str(num_files) ')']);

    try

        % run the file and store results, capturing all printed output
        [test_info, test_pass] = evalc([fn '(' plot_simulations ',' plot_comparisons ');']);
        % print the captured output
        fprintf('%s', test_info);

    catch %#ok<CTCH>
       
        % if the test gives an error for any reason, assign as failed
        test_pass = false;
        
    end
    
    % store test result
    test_result(filename_index) = test_pass;
    
end

% =========================================================================
% CREATE OUTPUT
% =========================================================================

completion_time = scaleTime(etime(clock, regression_start_time));
comp_info = getComputerInfo;
info = comp_info;
info.completion_time = completion_time;

% create results struct
test_struct = struct( ...
    'info', info, ...
    'results', struct('test', filenames(:), 'pass', num2cell(test_result(:)), 'test_info', test_info(:)) ...
);