% Replace broken MATLAB editor/run HTML links with plain text instructions.
%
% To address issue #5
% 
% In short, MATLAB 2024b+ opens help files in the system browser, breaking direct links.
% This script updates HTML files to show manual command instructions instead.

rootDir = pwd; % or specify your root directory
htmlFiles = dir(fullfile(rootDir, '**', '*.html'));

for k = 1:numel(htmlFiles)
    filePath = fullfile(htmlFiles(k).folder, htmlFiles(k).name);
    txt = fileread(filePath);

    % Replace edit links
    editPattern = '<li><a href="matlab:edit\(\[getkWavePath\(''([^'']*)''\) ''([^'']+?\.m)''\]\);" target="_top">Open the file in the MATLAB Editor</a></li>';
    editReplace = sprintf(['To open the file in MATLAB, enter the following command in the MATLAB Command Window:\n\nedit(getkWavePath(''$1'') ''$2'');\n']);

    txt = regexprep(txt, editPattern, editReplace, 'dotall');

    % Replace run links (with or without .m extension)
    runPattern = '<li><a href="matlab:run\(\[getkWavePath\(''([^'']*)''\) ''([^'']+?)(\.m)?''\]\);" target="_top">Run the file in MATLAB</a></li>';
    runReplace = sprintf(['To run the file, enter the following command in the MATLAB Command Window:\n\nrun(getkWavePath(''$1'') ''$2'');\n']);

    txt = regexprep(txt, runPattern, runReplace, 'dotall');

    % Only write if changes were made
    if ~strcmp(fileread(filePath), txt)
        fid = fopen(filePath, 'w');
        if fid == -1
            warning('Could not open file for writing: %s', filePath);
        else
            fwrite(fid, txt);
            fclose(fid);
            fprintf('Updated: %s\n', filePath);
        end
    end
end