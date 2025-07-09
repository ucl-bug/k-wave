function runUnitTests_show_results(test_struct)
%RUNUNITTESTS_ACTIONS Run MATLAB unit tests and format results for GitHub Actions.
%
% DESCRIPTION:
%     runUnitTests_actions processes the provided test_struct, saves the results
%     as a test_results.json artifact.
%

% =========================================================================
% DISPLAY SUMMARY
% =========================================================================

info = test_struct.info

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
            disp(fn);
            if ~isempty(test_info)
                fprintf(test_info);
                disp('  ');
            end
        end
    end
end