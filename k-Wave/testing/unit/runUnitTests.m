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

% remove this file from the list
filenames(contains(filenames, 'runUnitTests.m')) = [];

% filter filenames based on wildcard
if nargin > 0
    filenames = filenames(contains(filenames, wildcard));
end

% test just two tests
filenames = filenames(1:2)

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
% DISPLAY SUMMARY
% =========================================================================

% get information about PC
comp_info = getComputerInfo;
completion_time = scaleTime(etime(clock, regression_start_time));

% get k-Wave version
eval('cur_dir = pwd; cd(getkWavePath(''private'')); kwave_ver = getkWaveVersion; cd(cur_dir);');

% display test header
disp('   ');
disp('-------------------------------------------------------------------------------------');
disp('            _      __        __                _____         _            ');
disp('           | | __  \ \      / /_ ___   _____  |_   _|__  ___| |_ ___ _ __ ');
disp('           | |/ /___\ \ /\ / / _` \ \ / / _ \   | |/ _ \/ __| __/ _ \ ''__|');
disp('           |   <_____\ V  V / (_| |\ V /  __/   | |  __/\__ \ ||  __/ |   ');
disp('           |_|\_\     \_/\_/ \__,_| \_/ \___|   |_|\___||___/\__\___|_|   ');
disp('  ');                                                                
disp('-------------------------------------------------------------------------------------');
disp('  ');
disp(['DATE:                     ' comp_info.date]);
disp(['HOST NAME:                ' comp_info.computer_name]);
disp(['USER NAME:                ' comp_info.user_name]);
disp(['O/S TYPE:                 ' comp_info.operating_system_type]);
disp(['O/S:                      ' comp_info.operating_system]);
disp(['MATLAB VERSION:           ' comp_info.matlab_version]);
disp(['TESTED K-WAVE VERSION:    ' comp_info.kwave_version]);
disp(['TESTS COMPLETED IN:       ' completion_time]);
disp('  ');

% display individual test results
disp('UNIT TEST RESULTS:');

for filename_index = 1:length(filenames)
    
    % trim the filename
    fn = filenames{filename_index};
    fn = [fn(1:end - 2), ':'];
    
    % add some spaces to align results
    fn = sprintf('%-70s', fn);
    
    % append the test result
    if test_result(filename_index)
        disp(['  ' fn 'passed']);
    else
        disp(['  ' fn 'failed']);
    end
    
end

% display test summary
disp('  ');
disp('UNIT TEST SUMMARY:');
disp(['✅ Number of tests passed: ' num2str(sum(test_result))]);
disp(['❌ Number of tests failed: ' num2str( numel(test_result) - sum(test_result))]);
disp('  ');

if any(~test_result)
    disp('❌ FAILED TESTS:');
    for filename_index = 1:length(filenames)
        if ~test_result(filename_index)
            fn = filenames{filename_index};
            fn = fn(1:end - 2);
            disp(['--- ' fn ' ---']);
            fprintf('%s\n', test_info);
        end
    end
end

% =========================================================================
% CREATE OUTPUT
% =========================================================================

% add completion_time to comp_info
comp_info.completion_time = completion_time;

% create results struct
test_struct = struct( ...
    'info', comp_info, ...
    'results', struct('test', filenames(:), 'pass', num2cell(test_result(:)), 'test_info', test_info(:)) ...
);