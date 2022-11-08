function saveTiffStack(M, filename, bit_depth)
%SAVETIFFSTACK Save volume data as a tiff stack.
%
% DESCRIPTION:
%     saveTiffStack saves 3D volume data in matrix form as a stacked tiff
%     image. If the data is indexed as [x, y, z], each image corresponds to
%     an [y, z] slice of the data. If no filename is given, a dialog box is
%     invoked.
%
% USAGE:
%     saveTiffStack(M)
%     saveTiffStack(M, filename)
%     saveTiffStack(M, filename, bit_depth)
%     saveTiffStack(M, [], bit_depth)
%
% INPUTS:
%     M           - the volume data to save
%
% OPTIONAL INPUTS:
%     filename    - filename string
%     bit_depth   - bit depth of the tiff image (set to 8 or 16; default = 8) 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 15th June 2010
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby
%
% See also imwrite

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

% default bit depth
BIT_DEPTH_DEF = 8;

% set the bit depth of the image
if nargin < 3
    bit_depth = BIT_DEPTH_DEF;
end

% get filename if one isn't provided
if nargin == 1 || isempty(filename)
    filename = ['volumedata-' getDateString '.tif'];
    [file, path] = uiputfile(filename, 'Save file name');
    if isempty(file)
        error('File not selected.');
    end
    filename = [path, file];
end
    
% fail gracefully if no filename is selected
if ~isempty(filename)

    % update the command line status
    disp('Saving volume data as tiff stack...');

    % start the timer
    tic;

    % scale image data to maximise dynamic range
    scale_min = 0;
    scale_max = 2^bit_depth - 1;
    M = ( (M - min(M(:))) ./ (max(M(:) - min(M(:)))) ) * (scale_max - scale_min) + scale_min;
    switch bit_depth 
        case 8
            M = uint8(M);
        case 16
            M = uint16(M);
        otherwise
            error('Invalid bit depth.');
    end

    % write the last image frame
    imwrite(squeeze(M(end, :, :)), filename, 'tif', 'WriteMode', 'overwrite');

    % write the rest of the stack in reverse order
    for i=2:length(M(:, 1, 1))
        imwrite(squeeze(M(end - i + 1, :, :)), filename, 'tif', 'WriteMode', 'append');
    end

    % update command line status
    disp(['  completed in ' scaleTime(toc)]); 
end