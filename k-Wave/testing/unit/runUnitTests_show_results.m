function runUnitTests_show_results(test_struct)
%RUNUNITTESTS_SHOW_RESULTS Display MATLAB unit test results in a formatted summary.
%
% DESCRIPTION:
%     runUnitTests_show_results displays the results from the provided test_struct.
%
% =========================================================================
% DISPLAY SUMMARY
% =========================================================================

info = test_struct.info;
results = test_struct.results;

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
disp(['DATE:                     ' info.date]);
disp(['HOST NAME:                ' info.computer_name]);
disp(['USER NAME:                ' info.user_name]);
disp(['O/S TYPE:                 ' info.operating_system_type]);
disp(['O/S:                      ' info.operating_system]);
disp(['MATLAB VERSION:           ' info.matlab_version]);
disp(['TESTED K-WAVE VERSION:    ' info.kwave_version]);
disp(['TESTS COMPLETED IN:       ' info.completion_time]);
disp('  ');

% display test summary
disp('  ');
num_passed = sum([results.pass]);
num_failed = numel(results) - num_passed;
disp(['✅ Number of tests passed: ' num2str(num_passed)]);
disp(['❌ Number of tests failed: ' num2str(num_failed)]);
disp('  ');

% Show failed tests using test_struct
failed_idx = find(~[results.pass]);
if ~isempty(failed_idx)
    disp('FAILED TESTS:');
    for i = failed_idx
        fn = results(i).test;
        fn = fn(1:end - 2);
        disp(['❌  ' fn 'failed']);
    end
end

% display individual test results
disp('UNIT TEST RESULTS:');

for i = 1:length(results)
    
    % trim the filename
    fn = results(i).test;
    fn = [fn(1:end - 2), ':'];
    
    % add some spaces to align results
    fn = sprintf('%-70s', fn);
    
    % append the test result
    if results(i).pass
        disp(['✅  ' fn 'passed']);
    else
        disp(['❌  ' fn 'failed']);
    end
    
end

% display test summary
disp('NOTE:');
disp('Test output details are in the "Run unit tests" section of the workflow.');
disp('You can also download a JSON summary from the "Upload Artifact" section in your CI logs or dashboard.');

% Fail if any tests failed (for CI integration)
if num_failed > 0
    exit(1);
end
