function number_differences = h5compare(filename1, filename2)
%H5COMPARE Compare the contents of two HDF5 files
%
% DESCRIPTION: 
%     h5compare compares the datasets and attributes within two HDF5 files.
%     Fields that are present in one file but not the other are reported as
%     missing. For fields that exist in both files but have different
%     values, the L_inf error is reported. All outputs are printed to the
%     command line.
%
%     This function requires MATLAB 2011a or later.
%
% USAGE:
%     h5compare(filename1, filename2)
%
% INPUTS:
%     filename1             - filename of the first HDF5 file
%     filename2             - filaname of the second HDF5 file
%
% OUTPUTS:
%     number_differences    - number of differences between the files
%
% ABOUT:
%     author                - Bradley Treeby
%     date                  - 27th August 2012
%     last update           - 7th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby
%
% See also h5write, writeMatrix

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

% get the contents of files
info1 = h5info(filename1);
info2 = h5info(filename2);

% get the variables names
var_names1 = {info1.Datasets.Name};
var_names2 = {info2.Datasets.Name};

% create a combined list
var_names = unique([var_names1, var_names2]);

% set ok flag
number_differences = 0;

disp('---------------------------------------------------------');
disp('Comparing Datasets');
disp('---------------------------------------------------------'); 

% loop through each of the variables, then check the contents
for var_index = 1:length(var_names)
    
    % try to load data from file 1
    try
        
        % load file 1
        var1 = double(h5read(filename1, ['/' var_names{var_index}]));
        
        % try to load data from file 2
        try
            
            % load file 2
            var2 = double(h5read(filename2, ['/' var_names{var_index}]));

            % check values
            if ( ndims(var1) ~= ndims(var2) ) || ( any(size(var1) ~= size(var2)) )
                disp(['WARNING: Sizes for ' var_names{var_index} ' are different']);
                number_differences = number_differences + 1;
            elseif sum(abs(var1(:) - var2(:)))
                mx_error = max(abs(var1(:) - var2(:)));
                if  mx_error < 10*eps
                    disp(['Contents for ' var_names{var_index} ' are within 10^-15']);
                else
                    disp(['WARNING: Contents for ' var_names{var_index} ' are different (L_inf = ' num2str(mx_error) ')']);
                    number_differences = number_differences + 1;
                end
            else
                disp(['Contents for ' var_names{var_index} ' are equal']);
            end
            
        catch ME
            if strcmp(ME.identifier, 'MATLAB:imagesci:h5read:datasetDoesNotExist')  || strcmp(ME.identifier, 'MATLAB:imagesci:h5read:libraryError')
                disp(['WARNING: Contents for ' var_names{var_index} ' is missing in file 2']);
                number_differences = number_differences + 1;
            else
                rethrow(ME);
            end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5read:datasetDoesNotExist') || strcmp(ME.identifier, 'MATLAB:imagesci:h5read:libraryError')
            disp(['WARNING: Contents for ' var_names{var_index} ' is missing in file 1']);
            number_differences = number_differences + 1;
        else            
            rethrow(ME);
        end
    end
end

disp('---------------------------------------------------------');
disp('Comparing Attributes');
disp('---------------------------------------------------------'); 

% get the variables names
att_names = {info1.Attributes.Name};

% loop through each of the attributes, then check the contents
for att_index = 1:length(att_names)
    
    % load reference data
    att1 = h5readatt(filename1, '/', att_names{att_index});
    
    try
        
        % load comparison data
        att2 = h5readatt(filename2, '/', att_names{att_index});
        
        % check values
        if ~strcmp(att1, att2)
            disp(['WARNING: Contents for ' att_names{att_index} ' are different']);
            number_differences = number_differences + 1;
        else
            disp(['Contents for ' att_names{att_index} ' are equal']);
        end
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5readatt:cannotOpenAttribute')
            disp(['WARNING: Contents for ' att_names{att_index} ' is missing in second file']);
            number_differences = number_differences + 1;
        else
            rethrow(ME);
        end
    end
end

% display summary
disp('---------------------------------------------------------');
disp(['Total of ' num2str(number_differences) ' differences between the files']);
disp('---------------------------------------------------------');    