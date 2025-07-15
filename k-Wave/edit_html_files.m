
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

    % Pattern to match a <p> containing a <ul> with edit/run links (case-insensitive, multiline)
    pattern = ['<p>\s*<ul>\s*' ...
        '<li><a href="matlab:edit\(\[getkWavePath\(''([^'']*)''\) ''([^'']+?\.m)''\]\);" target="_top">Open the file in the MATLAB Editor</a></li>\s*' ...
        '<li><a href="matlab:run\(\[getkWavePath\(''([^'']*)''\) ''([^'']+?)''\]\);" target="_top">Run the file in MATLAB</a></li>\s*' ...
        '</ul>\s*</p>'];

    % Replacement with <p> preserved, but <ul>, <li>, and <a> replaced with instructions and code blocks
    replace = [newline ...
        'To open the file in MATLAB, enter the following command in the MATLAB Command Window.' newline ...
        '<pre class="codeinput">' newline ...
        'edit(getkWavePath(''$1'') ''$2'');' newline ...
        '</pre>' newline ...
        'To run the file, enter the following command in the MATLAB Command Window.' newline ...
        '<pre class="codeinput">' newline ...
        'run(getkWavePath(''$3'') ''$4'');' newline ...
        '</pre>' newline ...
        '</p>'];

    % Do the replacement
    newtxt = regexprep(txt, pattern, replace, 'dotall');

    % Only write if changes were made
    if ~strcmp(txt, newtxt)
        fid = fopen(filePath, 'w');
        if fid == -1
            warning('Could not open file for writing: %s', filePath);
        else
            fwrite(fid, newtxt);
            fclose(fid);
            fprintf('Updated: %s\n', filePath);
        end
    end
end