% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to set the
%     reference sound speed.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 25th June 2015
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2018 Bradley Treeby

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

% select the reference sound speed used in the k-space operator based on
% the heterogeneous sound speed map
if ~flags.elastic_code
    
    % calculate the reference sound speed for the fluid code, using the
    % maximum by default which ensures the model is unconditionally stable
    if isfield(medium, 'sound_speed_ref')
        if isnumeric(medium.sound_speed_ref)
            c_ref = medium.sound_speed_ref;
        elseif strcmp(medium.sound_speed_ref, 'min')
            c_ref = min(medium.sound_speed(:));
        elseif strcmp(medium.sound_speed_ref, 'mean')
            c_ref = mean(medium.sound_speed(:));
        elseif strcmp(medium.sound_speed_ref, 'max')
            c_ref = max(medium.sound_speed(:));        
        else
            error('Unknown input for medium.sound_speed_ref.');
        end
    else
        c_ref = max(medium.sound_speed(:));
    end
    disp(['  reference sound speed: ' num2str(c_ref) 'm/s']);
    
elseif ~flags.kspace_elastic_code

    % in the pstd elastic case, the reference sound speed is only used to
    % calculate the PML absorption, so just use the compressional wave
    % speed 
    if isfield(medium, 'sound_speed_ref')
        if isnumeric(medium.sound_speed_ref)
            c_ref = medium.sound_speed_ref;
        elseif strcmp(medium.sound_speed_ref, 'min')
            c_ref = min(medium.sound_speed_compression(:));
        elseif strcmp(medium.sound_speed_ref, 'mean')
            c_ref = mean(medium.sound_speed_compression(:));
        elseif strcmp(medium.sound_speed_ref, 'max')
            c_ref = max(medium.sound_speed_compression(:));
        else
            error('Unknown input for medium.sound_speed_ref.');
        end
    else
        c_ref = max(medium.sound_speed_compression(:));
    end
    disp(['  reference sound speed: ' num2str(c_ref) 'm/s']);    
    
else
    
    % in the k-space elastic case, there are two reference sound speeds for
    % the compressional and shear waves, so compute them seperately
    if isfield(medium, 'sound_speed_ref_compression')
        if isnumeric(medium.sound_speed_ref_compression)
            c_ref_compression = medium.sound_speed_ref_compression;
        elseif strcmp(medium.sound_speed_ref_compression, 'min')
            c_ref_compression = min(medium.sound_speed_compression(:));
        elseif strcmp(medium.sound_speed_ref_compression, 'mean')
            c_ref_compression = mean(medium.sound_speed_compression(:));
        elseif strcmp(medium.sound_speed_ref_compression, 'max')
            c_ref_compression = max(medium.sound_speed_compression(:));
        else
            error('Unknown input for medium.sound_speed_ref_compression.');
        end
    else
        c_ref_compression = max(medium.sound_speed_compression(:));
    end
    disp(['  reference sound speed (compression): ' num2str(c_ref_compression) 'm/s']);
    
    if isfield(medium, 'sound_speed_ref_shear')
        if isnumeric(medium.sound_speed_ref_shear)
            c_ref_shear = medium.sound_speed_ref_shear;
        elseif strcmp(medium.sound_speed_ref_shear, 'min')
            c_ref_shear = min(medium.sound_speed_shear(:));
        elseif strcmp(medium.sound_speed_ref_shear, 'mean')
            c_ref_shear = mean(medium.sound_speed_shear(:));
        elseif strcmp(medium.sound_speed_ref_shear, 'max')
            c_ref_shear = max(medium.sound_speed_shear(:));
        else
            error('Unknown input for medium.sound_speed_ref_shear.');
        end
    else
        c_ref_shear = max(medium.sound_speed_shear(:));
    end
    disp(['  reference sound speed (shear): ' num2str(c_ref_shear) 'm/s']);    
    
end