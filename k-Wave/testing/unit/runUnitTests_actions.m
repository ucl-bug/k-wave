function runUnitTests_actions(test_struct)
%RUNUNITTESTS_ACTIONS Run MATLAB unit tests and format results for GitHub Actions.
%
% DESCRIPTION:
%     runUnitTests_actions processes the provided test_struct, saves the results
%     as a test_results.json artifact.
%

% =========================================================================
% CREATE ARTIFACT
% =========================================================================

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
