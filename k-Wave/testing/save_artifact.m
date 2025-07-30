function save_artifact(test_struct)
%SAVE_ARTIFACT 
%
% DESCRIPTION:
%     save_artifact processes the provided test_struct, saves the results
%     as a test_results.json artifact.
%

fid = fopen('test_results.json', 'w');
if fid == -1
    warning('Could not open test_results.json for writing.');
else
    fwrite(fid, jsonencode(test_struct), 'char');
    fclose(fid);
end

