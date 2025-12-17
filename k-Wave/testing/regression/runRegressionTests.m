function test_struct = runRegressionTests(data_folder, show_results)
%RUNREGRESSIONTESTS Run MATLAB regression tests.
%
% DESCRIPTION:
%     runRegressionTests runs regression tests for the .mat files stored
%     within the same folder as runRegressionTests, or in the folder
%     specified by data_folder. The function works by generating a list of
%     the .mat files, then running the .m scripts with the same name (these
%     are assumed to be on the MATLAB path, or within the test_folder which
%     can be set below). The variables within the .mat file are then
%     compared to the variables in the workspace after the script has
%     executed. If the variable is a structure, the individual fields are
%     compared.
%
% USAGE:
%     runRegressionTests()
%     runRegressionTests(data_folder)
%
% INPUTS:
%     data_folder   - string of folder containing regression data
%     show_results  - logical flag to display results (default: true)
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 17th February 2014
%     last update   - 10th April 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014- Bradley Treeby

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

%#ok<*NODEF>
%#ok<*NASGU>
%#ok<*IDISVAR>

% Set defaults for optional arguments
if nargin < 2 || isempty(show_results)
    show_results = true;
end

% start the timer
regression_start_time = datetime('now');

% check for data_folder input and change directories
if nargin == 1
    cd(data_folder);
end

% set pass threshold
DOUBLE_COMPARISON_THRESHOLD = 5e-13;
SINGLE_COMPARISON_THRESHOLD = 1e-5;

% set folder for the .m files to be tested (leave blank if the .m files are
% on the MATLAB path)
regression_test_folder = getkWavePath('examples');

% clear previous regression test results (this filename is used explicitly
% later in the file to circumvent clear all)
temp_var_filename = [tempdir, 'runRegressionTests_TEMP_VARS.mat'];
if exist(temp_var_filename, 'file')
    delete(temp_var_filename);
end

% get a list of the .mat files
filenames = what;
filenames = filenames.mat;

% extract number of regression tests to run
num_files = length(filenames);

% keep a list of whether the test passed or failed
test_result = false(num_files, 1);

% preallocate cell array for test_info with empty strings by default
test_info = repmat({''}, num_files, 1);

% =========================================================================
% RUN TESTS
% =========================================================================

% run the examples one by one
for filename_index = 1:num_files

    % trim the .mat extension to run the corresponding m-file
    fn = filenames{filename_index};
    fn = fn(1:end-4);
    
    % save required variables to the workspace to circumvent clear all
    % within the example files
    save(temp_var_filename, ...
        'DOUBLE_COMPARISON_THRESHOLD', 'SINGLE_COMPARISON_THRESHOLD', ...
        'regression_test_folder', 'temp_var_filename', 'regression_start_time', ...
        'filename_index', 'filenames', 'fn', 'num_files', 'test_result', 'show_results');
    
    % display the filename
    disp(' ');
    disp(['Running regression test for ' fn ' (Test ' num2str(filename_index) ' of ' num2str(num_files) ')']);

    % run the example file
    run([regression_test_folder, fn]);
    
    % close the figure windows
    close all;

    % load back the workspace variables cleared by the example file
    load([tempdir, 'runRegressionTests_TEMP_VARS.mat']); %#ok<LOAD>
    
    % get the list of variables in the corresponding regression data
    vars = whos('-file', fn);
        
    % assign comparison threshold
    load(fn, 'comp_info');
    if strcmp(comp_info.precision, 'single')
        comparison_thresh = SINGLE_COMPARISON_THRESHOLD;
    elseif strcmp(comp_info.precision, 'double')
        comparison_thresh = DOUBLE_COMPARISON_THRESHOLD;
    else
        error('unknown precision setting');
    end
    
    % assign test pass variable
    test_pass_overall = true;
    
    % loop through the variables
    for var_index = 1:length(vars)
        
        % check variable is not the comp_info variable
        if ~strcmp(vars(var_index).name, 'comp_info')
            
            % clear some things
            clear reg_data;
            
            % load the variable to the workspace
            reg_data = load(fn, vars(var_index).name); 
            
            % check if variable is a structure or an object
            if isstruct(eval(['reg_data.' vars(var_index).name])) || isobject(eval(['reg_data.' vars(var_index).name]))

                % get the field names
                field_names = fieldnames(eval(['reg_data.' vars(var_index).name]));

                % loop through each of the structure fields
                for field_index = 1:length(field_names)
                   
                    % display what's happening
                    fprintf(['comparing ' vars(var_index).name '.' field_names{field_index} '... ']);
                    
                    % load the regression data
                    eval(['var_ref = reg_data.' vars(var_index).name '.' field_names{field_index} ';']);

                    % give a pseudonym to the current data
                    eval(['var_new = ' vars(var_index).name '.' field_names{field_index} ';']);

                    % check size, and then compare error
                    if any(size(var_new) ~= size(var_ref))
                        test_pass = false;
                        fprintf('failed (data sizes different)\n');
                    else
                    
                        % only compare if not empty, or inf
                        if isempty(var_ref) && isempty(var_new)
                            test_pass = true;
                            fprintf('passed (empty)\n');
                        elseif isempty(var_ref) || isempty(var_new)
                            test_pass = false;
                            fprintf('failed (one of the variables empty)\n');
                        elseif (numel(var_ref) == 1) && (numel(var_new) == 1) && isinf(var_ref) && isinf(var_new)
                            test_pass = true;
                            fprintf('passed (inf)\n');
                        elseif (numel(var_ref) == 1) && (numel(var_new) == 1) && (isinf(var_ref) || isinf(var_new))
                            test_pass = false;
                            fprintf('failed (one of the variables inf)\n');                            
                        else
                        
                            % compute error metric, avoiding divide by 0
                            if max(abs(var_ref(:))) == 0 
                                L_inf = max(abs(var_new(:) - var_ref(:)));
                            else
                                L_inf = max(abs(var_new(:) - var_ref(:))) / max(abs(var_ref(:)));
                            end

                            % check if test passed
                            if L_inf > comparison_thresh || isnan(L_inf)
                                test_pass = false;
                                fprintf('failed');
                            else
                                test_pass = true;
                                fprintf('passed');
                            end

                            % display the error
                            fprintf(' (L_inf = %e)\n', L_inf);
                            
                        end
                        
                    end
                    
                    % clear the variables just in case
                    clear var_ref var_new;
                    
                    % store overall test pass
                    test_pass_overall = test_pass_overall && test_pass;

                end

            else

                % display what's happening
                fprintf(['comparing ' vars(var_index).name '... ']);
                
                % load the regression data
                eval(['var_ref = reg_data.' vars(var_index).name ';']);

                % give a pseudonym to the current data
                eval(['var_new = ' vars(var_index).name ';']);

                % check size, and then compare error
                if any(size(var_new) ~= size(var_ref))
                    test_pass = false;
                    fprintf('failed (data sizes different)\n');
                else                

                    % compute error metric, avoiding divide by 0
                    if max(abs(var_ref(:))) == 0
                        L_inf = max(abs(var_new(:) - var_ref(:)));
                    else
                        L_inf = max(abs(var_new(:) - var_ref(:))) / max(abs(var_ref(:)));
                    end

                    % check if test passed
                    if L_inf > comparison_thresh || isnan(L_inf)
                        test_pass = false;
                        fprintf('failed');
                    else
                        test_pass = true;
                        fprintf('passed');
                    end
                    % save the error string
                    error_str = sprintf(' (L_inf = %e)\n', L_inf);

                    % display the error
                    fprintf('%s', error_str);

                    % save L_inf as test_info
                    test_info{filename_index} = error_str;
                end
                
                % clear the variables just in case
                clear var_ref var_new;
                
                % store overall test pass
                test_pass_overall = test_pass_overall && test_pass;

            end
        end
    end
    
    % store test result
    test_result(filename_index) = test_pass_overall;
end

% =========================================================================
% CREATE OUTPUT
% =========================================================================

completion_time = scaleTime(seconds(datetime('now') - regression_start_time));
comp_info = getComputerInfo;
info = comp_info;
info.completion_time = completion_time;

% create results struct
test_struct = struct( ...
    'info', info, ...
    'results', struct('test', filenames(:), 'pass', num2cell(test_result(:)), 'test_info', test_info(:)) ...
);

% =========================================================================
% SHOW RESULTS
% =========================================================================

if show_results
    show_test_results(test_struct);
end
