function  runUnitTests_actions(test_struct)
%RUNUNITTESTS_ACTIONS Run MATLAB unit tests and format results for GitHub Actions.
%
% DESCRIPTION:
%     runUnitTests_actions runs all unit tests and processes the results for
%     GitHub Actions. It creates a test_results.json artifact and displays
%     the results as a markdown table for easy viewing in the GitHub Actions log.
%

% =========================================================================
% DISPLAY AS MD TABLE
% =========================================================================

n_tests = size(test_struct.results,1);
n_passed = sum([test_struct.results.pass]);
n_failed = n_tests - n_passed;

disp('**Test Environment Info**');
info_fields = fieldnames(test_struct.info);
for i = 1:numel(info_fields)
    field = info_fields{i};
    value = test_struct.info.(field);
    if isnumeric(value)
        value = num2str(value);
    end
    fprintf('%s: %s\n', field, value);
end
disp(' ');
disp('**Test Results Table**');
disp('> Click on a test name to expand and view its output or error details.');
disp(' ');
disp('|   | Test Name |');
disp('|---|-----------|');

for i = 1:n_tests
    if test_struct.results(i).pass
        status = '✅';
    else
        status = '❌';
    end
    % Test name with foldout menu for details
    if isfield(test_struct.results, 'test_info') && numel(test_struct.results) >= i
        info = test_struct.results(i).test_info;
    else
        info = '';
    end
    test_name = test_struct.results(i).test;
    test_name_md = sprintf('<details><summary>%s</summary>\n\n```\n%s\n```\n</details>', test_name, info);
    % Print Markdown row
    fprintf('| %s | %s |\n', status, test_name_md);
end

disp(' ');
disp('---');
disp('**Test Summary**');
fprintf('✅ Passed: %d\n', n_passed);
fprintf('❌ Failed: %d\n', n_failed);
disp(' ');


% =========================================================================
% CREATE ARTIFACT
% =========================================================================

% Save to json file
fid = fopen('test_results.json', 'w');
if fid == -1
    warning('Could not open test_results.json for writing.');
else
    fwrite(fid, jsonencode(test_struct), 'char');
    fclose(fid);
end

disp(' ');
disp('UNIT TEST RESULTS SAVED TO test_results.json');
disp(' ');
