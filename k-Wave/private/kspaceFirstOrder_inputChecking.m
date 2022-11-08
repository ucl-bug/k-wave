% DESCRIPTION:
%     Subscript to check input structures and optional input parameters.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 21st December 2010
%     last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby

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

%#ok<*NASGU>

% =========================================================================
% CHECK DIMENSIONS AND FUNCTION
% =========================================================================

% get the name of the calling function
calling_func = dbstack;

% check for deprecated field_data outputs
if num_outputs == 2 && ~strcmp(calling_func(3).name, 'benchmark')
    error(['The output field_data has been deprecated. Please use sensor.record = {''p_final''} (etc) to return the final pressure field. See ' calling_func(2).name ' documentation for more information.']);
end

% check correct function has been called for the dimensionality of kgrid
switch calling_func(2).name
    case 'kspaceFirstOrder1D'
        if kgrid.dim ~= 1
            error(['kgrid has the wrong dimensionality for ' calling_func(2).name '.']);
        end
    case {'kspaceFirstOrder2D', 'pstdElastic2D', 'kspaceElastic2D', 'kspaceFirstOrderAS'}
        if kgrid.dim ~= 2
            error(['kgrid has the wrong dimensionality for ' calling_func(2).name '.']);
        end        
    case {'kspaceFirstOrder3D', 'pstdElastic3D', 'kspaceElastic3D'}
        if kgrid.dim ~= 3
            error(['kgrid has the wrong dimensionality for ' calling_func(2).name '.']);
        end        
end

% run subscript to set defaults for flags and optional inputs
kspaceFirstOrder_setDefaults;

% check whether the calling function is a fluid or elastic code
if strncmp(calling_func(2).name, 'pstdElastic', 11) || strncmp(calling_func(2).name, 'kspaceElastic', 13)
    
    % set elastic code flag
    flags.elastic_code = true;
    
    % set additional flag if elastic code is a k-space code, as this
    % requires some additional variables
    if strncmp(calling_func(2).name, 'kspaceElastic', 13)
        flags.kspace_elastic_code = true;
    end

end

% check whether the code is axisymmetric
if strcmp(calling_func(2).name, 'kspaceFirstOrderAS')
    flags.axisymmetric = true;
end

% cleanup unused variables
clear calling_func;

% update command line status with the start time
if flags.elastic_code
    disp('Running k-Wave elastic simulation...');
else
    disp('Running k-Wave simulation...');
end
disp(['  start time: ' datestr(start_time)]);

% =========================================================================
% INDEX DATA TYPES
% =========================================================================

% pre-calculate the data type needed to store the matrix indices given the
% total number of grid points: indexing variables will be created using
% this data type to save memory
if kgrid.total_grid_points < intmax('uint8')
    index_data_type = 'uint8';
elseif kgrid.total_grid_points < intmax('uint16')
    index_data_type = 'uint16';
elseif kgrid.total_grid_points < intmax('uint32')
    index_data_type = 'uint32';                
else
    index_data_type = 'double';
end   

% =========================================================================
% CHECK OPTIONAL INPUTS
% =========================================================================

% run subscript to check optional inputs
kspaceFirstOrder_checkOptionalInputs;

% =========================================================================
% CHECK MEDIUM STRUCTURE INPUTS
% =========================================================================

% check the medium input is defined as a structure
if ~isstruct(medium)
    error('medium must be defined as a MATLAB structure.');
end

% check medium fields (correctly setting the list of allowable field names
% avoids needing additional checks later on, for example, whether BonA has
% been set for the elastic code which is not allowed)
if flags.elastic_code && flags.kspace_elastic_code
    
    % check the allowable field names
    checkFieldNames(medium, {'sound_speed_shear', 'sound_speed_compression', 'density',...
        'sound_speed_ref_compression', 'sound_speed_ref_shear',...
        'alpha_coeff_shear', 'alpha_coeff_compression', 'alpha_power_shear', 'alpha_power_compression'});
    
    % force the sound speed and density to be defined
    enforceFields(medium, {'sound_speed_compression', 'sound_speed_shear', 'density'});
    
elseif flags.elastic_code
    
    % check the allowable field names
    checkFieldNames(medium, {'sound_speed_shear', 'sound_speed_compression', 'sound_speed_ref', 'density',...
        'alpha_coeff_shear', 'alpha_coeff_compression'});
    
    % force the sound speed and density to be defined
    enforceFields(medium, {'sound_speed_compression', 'sound_speed_shear', 'density'});    
    
else
    
    % check the allowable field names
    checkFieldNames(medium, {'sound_speed', 'sound_speed_ref', 'density',...
        'alpha_coeff', 'alpha_power', 'alpha_mode', 'alpha_filter', 'alpha_sign', 'BonA'});
    
    % force the sound speed to be defined
    enforceFields(medium, {'sound_speed'});
    
end

% if using the fluid code, allow the density field to be blank if the
% medium is homogeneous
if ~flags.elastic_code && ~isfield(medium, 'density') && numel(medium.sound_speed) == 1
    user_medium_density_input = false;
    medium.density = 1;
else
    enforceFields(medium, {'density'});
    user_medium_density_input = true;
end

% check medium absorption inputs for the fluid code
if isfield(medium, 'alpha_coeff') || isfield(medium, 'alpha_power')
    
    % set absorption flag
    flags.absorbing = true;
    
    % check absorption mode (only stokes absorption is supported in the
    % axisymmetric code)
    if flags.axisymmetric || (isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'stokes'))
         
        % set equation of state variable
        equation_of_state = 'stokes';
        
        % set stokes absorption flag
        flags.stokes = true;
        
        % enforce absorption coefficient
        enforceFields(medium, {'alpha_coeff'});
        
        % give warning if y is specified
        if isfield(medium, 'alpha_power') && ((numel(medium.alpha_power) ~= 1) || (medium.alpha_power ~= 2))
            disp('  WARNING: the axisymmetric code and stokes absorption assume alpha_power = 2, user value ignored.');
        end
        
        % overwrite y value
        medium.alpha_power = 2; 
        
        % don't allow medium.alpha_mode with the axisymmetric code
        if isfield(medium, 'alpha_mode') && (strcmp(medium.alpha_mode, 'no_absorption') || strcmp(medium.alpha_mode, 'no_dispersion'))
            error('Input option medium.alpha_mode is not supported with the axisymmetric code or medium.alpha_mode = ''stokes''.');
        end
        
        % don't allow alpha_filter with stokes absorption (no variables are
        % applied in k-space)
        if isfield(medium, 'alpha_filter')
            error('Input option medium.alpha_filter is not supported with the axisymmetric code or medium.alpha_mode = ''stokes''.');
        end
        
    else

        % set equation of state variable
        equation_of_state = 'absorbing';
        
        % enforce both absorption parameters
        enforceFields(medium, {'alpha_coeff', 'alpha_power'});        
        
        % check y is a scalar
        if numel(medium.alpha_power) ~= 1
            error('medium.alpha_power must be scalar.');
        end

        % check y is real and within 0 to 3
        if ~isreal(medium.alpha_coeff) || (medium.alpha_power >= 3 || medium.alpha_power < 0)
            error('medium.alpha_power must be a real number between 0 and 3.');
        end

        % display warning if y is close to 1 and the dispersion term has
        % not been set to zero
        if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_dispersion'))
            if medium.alpha_power == 1
                error('The power law dispersion term in the equation of state is not valid for medium.alpha_power = 1. This error can be avoided by choosing a power law exponent close to, but not exactly, 1. If modelling acoustic absorption for medium.alpha_power = 1 is important and modelling dispersion is not critical, this error can also be avoided by setting medium.alpha_mode to ''no_dispersion''.');
            end
        end
        
    end
    
    % check the absorption mode input is valid
    if isfield(medium, 'alpha_mode')
        if ~ischar(medium.alpha_mode) || (~strcmp(medium.alpha_mode, 'no_absorption') && ~strcmp(medium.alpha_mode, 'no_dispersion') && ~strcmp(medium.alpha_mode, 'stokes'))
            error('medium.alpha_mode must be set to ''no_absorption'', ''no_dispersion'', or ''stokes''.');
        end
    end    
    
    % check the absorption filter input is valid
    if isfield(medium, 'alpha_filter') && ~all(size(medium.alpha_filter) == size(kgrid.k))
        error('medium.alpha_filter must be the same size as the computational grid.');
    end
    
    % check the absorption sign input is valid
    if isfield(medium, 'alpha_sign') && (~isnumeric(medium.alpha_sign) || (numel(medium.alpha_sign) ~= 2))
        error('medium.alpha_sign must be given as a 2 element numerical array controlling absorption and dispersion, respectively.')
    end
    
    % check alpha_coeff is non-negative and real
    if ~isreal(medium.alpha_coeff) || any(medium.alpha_coeff(:) < 0)
        error('medium.alpha_coeff must be non-negative and real.');
    end
    
else
    
    % set equation of state variable
    equation_of_state = 'lossless';
    
end

% check medium absorption inputs for the elastic code
if flags.elastic_code
    if isfield(medium, 'alpha_coeff_compression') || isfield(medium, 'alpha_coeff_shear')

        % if one absorption parameter is given, enforce the other
        enforceFields(medium, {'alpha_coeff_compression', 'alpha_coeff_shear'});

        % set elastic absorption flags
        flags.kelvin_voigt_model = true;
        flags.absorbing = true;

    end
end

% check if BonA is given and then set the nonlinear flag
if isfield(medium, 'BonA')
    flags.nonlinear = true;
end

% select the reference sound speed used in the k-space operator
kspaceFirstOrder_setSoundSpeedRef;

% =========================================================================
% CHECK SENSOR STRUCTURE INPUTS
% =========================================================================

% check sensor fields
if ~isempty(sensor)

    % check the sensor input is valid
    if ~(isstruct(sensor) || isa(sensor, 'kWaveTransducer'))
        error('sensor must be defined as a MATLAB structure or an object of the kWaveTransducer class.');
    end
    
    % set sensor flag to true
    flags.use_sensor = true;
    
    % check if sensor is a transducer, otherwise check input fields
    if ~isa(sensor, 'kWaveTransducer')
        if kgrid.dim == 2

            % check field names, including the directivity inputs for the
            % regular 2D code, but not the axisymmetric code
            if flags.axisymmetric
                checkFieldNames(sensor, {'mask', 'time_reversal_boundary_data', ...
                    'frequency_response', 'record_mode', 'record', 'record_start_index'});                
            else
                checkFieldNames(sensor, {'mask', 'directivity_pattern', 'directivity_angle', 'directivity_size',...
                    'time_reversal_boundary_data', 'frequency_response', 'record_mode', 'record', 'record_start_index'});
            end
                
            % check for sensor directivity input and set flag
            if isfield(sensor, 'directivity_angle')
                
                % set flag
                flags.compute_directivity = true;
                
                % make sure the sensor mask is not blank
                enforceFields(sensor, {'mask'});
                
                % check sensor.directivity_pattern and sensor.mask have the same size
                if sum(size(sensor.directivity_angle) ~= size(sensor.mask))
                    error('sensor.directivity_angle and sensor.mask must be the same size.')
                end
                
                % check if directivity pattern input exists, otherwise
                % apply default
                if ~isfield(sensor,'directivity_pattern')
                    sensor.directivity_pattern = DIRECTIVITY_PATTERN_DEF;
                end
                
                % check if directivity size input exists, otherwise make it
                % a constant times kgrid.dx
                if ~isfield(sensor,'directivity_size')
                    sensor.directivity_size = DIRECTIVITY_SIZE_DEF * max([kgrid.dx, kgrid.dy]);
                end
                
                % find the unique directivity angles
                sensor.directivity_unique_angles = unique(sensor.directivity_angle(sensor.mask == 1));
                
                % assign the wavenumber vectors
                sensor.directivity_wavenumbers = [kgrid.ky(:)'; kgrid.kx(:)'];

            end

        else
            
            % check field names without directivity inputs (these are not supported in 1 or 3D)
            checkFieldNames(sensor, {'mask', 'time_reversal_boundary_data', ...
                'frequency_response', 'record_mode', 'record', 'record_start_index'});
            
        end
        
        % check for time reversal inputs and set flags
        if ~flags.elastic_code && isfield(sensor, 'time_reversal_boundary_data')
            flags.time_rev = true;
            flags.record_p = false;
        end
        
        % check for the deprecated 'record_mode' input and throw error
        if isfield(sensor, 'record_mode')
            error('sensor.record_mode input has been deprecated. Please use the syntax sensor.record = {''p_rms'', ''p_max'', ...} to set recorded fields.');
        end
        
        % check for sensor.record and set usage flags - if no flags are
        % given, the time history of the acoustic pressure is recorded by
        % default 
        if isfield(sensor, 'record')
            
            % check for time reversal data
            if flags.time_rev
                disp('WARNING: sensor.record is not used for time reversal reconstructions');
            end
            
            % check the input is a cell array
            if ~iscell(sensor.record)
                error('sensor.record must be given as a cell array, e.g., {''p'', ''u''}.');
            end
           
            % check the sensor record flags
            if flags.elastic_code
                
                % list of allowable record flags for the elastic code
                record.flags = {'p', 'p_max', 'p_min', 'p_rms', 'p_max_all', 'p_min_all', 'p_final', ...
                    'u', 'u_max', 'u_min', 'u_rms', 'u_max_all', 'u_min_all', 'u_final', ...
                    'u_non_staggered', 'u_split_field', ...
                    'I', 'I_avg'};
                
                % check the contents of the cell array are valid inputs
                for record_index = 1:length(sensor.record)
                    if ~ismember(sensor.record(record_index), record.flags)
                        error(['''' sensor.record{record_index} ''' is not a valid input for sensor.record.']);
                    end
                end
                
            else
                
                % list of allowable record flags for the fluid code
                record.flags = {'p', 'p_max', 'p_min', 'p_rms', 'p_max_all', 'p_min_all', 'p_final', ...
                            'u', 'u_max', 'u_min', 'u_rms', 'u_max_all', 'u_min_all', 'u_final', ...
                            'u_non_staggered', ...
                            'I', 'I_avg'};   
                
                % check the contents of the cell array are valid inputs
                for record_index = 1:length(sensor.record)
                    if ~ismember(sensor.record(record_index), record.flags)
                        error(['''' sensor.record{record_index} ''' is not a valid input for sensor.record.']);
                    end
                end
                
            end
            
            % set the usage flags depending on the cell array
            if ismember('p', sensor.record)
                flags.record_p = true;
            else
                
                % set flags.record_p to false if a user input for sensor.record
                % is given and 'p' is not set (default is true)
                flags.record_p = false;
                
            end
            if ismember('p_max', sensor.record)
                flags.record_p_max = true;
            end   
            if ismember('p_min', sensor.record)
                flags.record_p_min = true;
            end             
            if ismember('p_rms', sensor.record)
                flags.record_p_rms = true;
            end 
            if ismember('p_max_all', sensor.record)
                flags.record_p_max_all = true;
            end   
            if ismember('p_min_all', sensor.record)
                flags.record_p_min_all = true;
            end             
            if ismember('p_final', sensor.record)
                flags.record_p_final = true;
            end  
            if ismember('u', sensor.record)
                flags.record_u = true;
            end            
            if ismember('u_max', sensor.record)
                flags.record_u_max = true;
            end
            if ismember('u_min', sensor.record)
                flags.record_u_min = true;
            end            
            if ismember('u_rms', sensor.record)
                flags.record_u_rms = true;
            end
            if ismember('u_max_all', sensor.record)
                flags.record_u_max_all = true;
            end
            if ismember('u_min_all', sensor.record)
                flags.record_u_min_all = true;
            end              
            if ismember('u_final', sensor.record)
                flags.record_u_final = true;
            end
            if ismember('u_non_staggered', sensor.record)
                flags.record_u_non_staggered = true;
            end  
            if ismember('u_split_field', sensor.record)
                flags.record_u_split_field = true;
            end              
            if ismember('I', sensor.record)
                flags.record_I = true;
            end     
            if ismember('I_avg', sensor.record)
                flags.record_I_avg = true;
            end
            
        end
        
        % enforce the sensor.mask field unless just recording the max_all
        % and _final variables
        if flags.record_p || flags.record_p_max || flags.record_p_min || flags.record_p_rms || ...
                flags.record_u || flags.record_u_non_staggered || flags.record_u_split_field || flags.record_u_max ||...
                flags.record_u_min || flags.record_u_rms || flags.record_I || flags.record_I_avg
            enforceFields(sensor, {'mask'});
        elseif ~flags.time_rev
            flags.blank_sensor = true;
        end
        
        % check for a user input for record_start_index
        if ~isfield(sensor, 'record_start_index')
            
            % if no input is given, record the time series from the
            % beginning 
            sensor.record_start_index = 1;
            
        else
            
            % force the user index to be an integer (this is checked later
            % after first defining kgrid.t_array if it is not defined)
            sensor.record_start_index = round(sensor.record_start_index);
            
        end
        
        % check if sensor mask is a binary grid, a set of cuboid corners,
        % or a set of Cartesian interpolation points
        if ~flags.blank_sensor
            if (kgrid.dim == 3 && numDim(sensor.mask) == 3) || (kgrid.dim ~= 3 && all(size(sensor.mask) == size(kgrid.k)))

                % check the grid is binary
                if sum(sensor.mask(:)) ~= numel(sensor.mask) - sum(sensor.mask(:) == 0)
                    error('sensor.mask must be a binary grid (numeric values must be 0 or 1).');
                end

                % check the grid is not empty
                if sum(sensor.mask(:)) == 0
                    error('sensor.mask must be a binary grid with at least one element set to 1.');
                end

            elseif size(sensor.mask, 1) == (2 * kgrid.dim)
                
                % set cuboid_corners flag
                flags.cuboid_corners = true;
                
                % make sure the points are integers
                if ~all(mod(sensor.mask(:), 1) == 0)
                    error('sensor.mask cuboid corner indices must be integers.');
                end
                
                % store a copy of the cuboid corners
                record.cuboid_corners_list = sensor.mask;
                
                % check the list makes sense
                if any(any(sensor.mask(1 + kgrid.dim:end, :) - sensor.mask(1:kgrid.dim, :) < 0))
                    switch kgrid.dim
                        case 1
                            error('sensor.mask cuboid corners must be defined as [x1, x2; ...].'' where x2 => x1, etc.');
                        case 2
                            error('sensor.mask cuboid corners must be defined as [x1, y1, x2, y2; ...].'' where x2 => x1, etc.');
                        case 3
                            error('sensor.mask cuboid corners must be defined as [x1, y1, z1, x2, y2, z2; ...].'' where x2 => x1, etc.');
                    end
                end
                
                % check the list are within bounds
                if any(sensor.mask < 1)
                    error('sensor.mask cuboid corners must be within the grid.');
                else
                    switch kgrid.dim
                        case 1
                            if any(sensor.mask > kgrid.Nx)
                                error('sensor.mask cuboid corners must be within the grid.');
                            end
                        case 2
                            if any(any(sensor.mask([1, 3], :) > kgrid.Nx)) || any(any(sensor.mask([2, 4], :) > kgrid.Ny))
                                error('sensor.mask cuboid corners must be within the grid.');
                            end
                        case 3
                            if any(any(any(sensor.mask([1, 4], :) > kgrid.Nx))) || any(any(any(sensor.mask([2, 5], :) > kgrid.Ny))) || any(any(any(sensor.mask([3, 6], :) > kgrid.Nz)))
                                error('sensor.mask cuboid corners must be within the grid.');
                            end
                    end
                end
                
                % create a binary mask for display from the list of corners
                sensor.mask = false(size(kgrid.k));
                for cuboid_index = 1:size(record.cuboid_corners_list, 2)
                    switch kgrid.dim
                        case 1
                            sensor.mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(2, cuboid_index)) = 1;
                        case 2
                            sensor.mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(3, cuboid_index),...
                                        record.cuboid_corners_list(2, cuboid_index):record.cuboid_corners_list(4, cuboid_index)) = 1;                       
                        case 3
                            sensor.mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(4, cuboid_index),...
                                        record.cuboid_corners_list(2, cuboid_index):record.cuboid_corners_list(5, cuboid_index),...
                                        record.cuboid_corners_list(3, cuboid_index):record.cuboid_corners_list(6, cuboid_index)) = 1;
                    end
                end
                
            else

                % check the Cartesian sensor mask is the correct size 
                % (1 x N, 2 x N, 3 x N)
                if size(sensor.mask, 1) ~= kgrid.dim || numDim(sensor.mask) > 2
                    error(['Cartesian sensor.mask for a ' num2str(kgrid.dim) 'D simulation must be given as a ' num2str(kgrid.dim) ' by N array.']);
                end

                % set Cartesian mask flag (this is modified in
                % createStorageVariables if the interpolation setting is
                % set to nearest)
                flags.binary_sensor_mask = false;   

                % extract Cartesian data from sensor mask
                switch kgrid.dim
                    case 1
                        
                        % align sensor data as a column vector to be the
                        % same as kgrid.x_vec so that calls to interp1
                        % return data in the correct dimension 
                        sensor_x = reshape(sensor.mask, [], 1);

                        % add sensor_x to the record structure for use with
                        % the _extractSensorData subfunction
                        record.sensor_x = sensor_x;

                    case 2
                        sensor_x = sensor.mask(1, :);
                        sensor_y = sensor.mask(2, :);
                    case 3
                        sensor_x = sensor.mask(1, :);
                        sensor_y = sensor.mask(2, :);
                        sensor_z = sensor.mask(3, :);    
                end

                % compute an equivalent sensor mask using nearest neighbour
                % interpolation, if flags.time_rev = false and
                % cartesian_interp = 'linear' then this is only used for
                % display, if flags.time_rev = true or cartesian_interp =
                % 'nearest' this grid is used as the sensor.mask  
                [sensor.mask, order_index, reorder_index] = cart2grid(kgrid, sensor.mask, flags.axisymmetric);

                % if in time reversal mode, reorder the p0 input data in
                % the order of the binary sensor_mask  
                if flags.time_rev

                    % append the reordering data
                    new_col_pos = length(sensor.time_reversal_boundary_data(1, :)) + 1;
                    sensor.time_reversal_boundary_data(:, new_col_pos) = order_index;

                    % reorder p0 based on the order_index
                    sensor.time_reversal_boundary_data = sortrows(sensor.time_reversal_boundary_data, new_col_pos);

                    % remove the reordering data
                    sensor.time_reversal_boundary_data = sensor.time_reversal_boundary_data(:, 1:new_col_pos - 1);

                end
            end
        end
        
    else

        % if the sensor is a transducer, check that the simulation is in 3D
        if kgrid.dim ~= 3
            error('Transducer inputs are only compatible with 3D simulations.');
        end
        
        % check that the transducer is only being used in forward mode
        if isfield(sensor, 'time_reversal_boundary_data')
            error('Transducer inputs not supported with time reversal.');
        end
        
        % check that sensor.record is not set
        if isfield(sensor, 'record')
            error('sensor.record cannot be set when using kWaveTransducer objects as the sensor.');
        end
        
        % set transducer sensor flag
        flags.transducer_sensor = true;
        flags.record_p = false;
        
        % check to see if there is an elevation focus
        if ~isinf(sensor.elevation_focus_distance)
            
            % set flag
            flags.transducer_receive_elevation_focus = true;
            
            % get the elevation mask that is used to extract the correct values
            % from the sensor data buffer for averaging
            transducer_receive_mask = sensor.elevation_beamforming_mask;
            
        end
        
    end
end

% if using time reversal with the elastic code, reassign the time reversal
% data as a pressure source (normal stress) with a dirichlet boundary
% condition. Note: this replaces the user inputs for source and sensor
if flags.elastic_code && isfield(sensor, 'time_reversal_boundary_data')
    
    % define a new source structure
    clear source;
    source.s_mask = sensor.mask;
    source.sxx = -flipdim(sensor.time_reversal_boundary_data, 2);
    if kgrid.dim >= 2
        source.syy = source.sxx;
    end
    if kgrid.dim == 3
        source.szz = source.sxx;
    end
    source.s_mode = 'dirichlet';
    
    % define a new sensor structure
    clear sensor;
    sensor.record = {'p_final'};
    
    % set flags for elastic time reversal
    flags.elastic_time_rev     = true;
    flags.use_sensor           = true;
    flags.blank_sensor         = true;
    sensor.record_start_index  = 1;
    
end

% check for directivity inputs with time reversal
if kgrid.dim == 2 && flags.use_sensor && flags.compute_directivity && flags.time_rev
    disp('WARNING: sensor directivity fields are not used for time reversal.');
end

% =========================================================================
% CHECK SOURCE STRUCTURE INPUTS
% =========================================================================

% check source inputs
if ~(isstruct(source) ||isa(source, 'kWaveTransducer'))
    
    % allow an invalid or empty source input if computing time reversal,
    % otherwise return error
    if ~flags.time_rev
        error('source must be defined as a MATLAB structure or an object of the kWaveTransducer class.');
    end
    
elseif ~isa(source, 'kWaveTransducer')
    
    % --------------------------
    % SOURCE IS NOT A TRANSDUCER
    % --------------------------
    
    % check allowable source types
    switch kgrid.dim
        case 1
            if ~flags.elastic_code
                checkFieldNames(source, {'p0', 'p', 'p_mask', 'p_mode', 'p_frequency_ref', 'ux', 'u_mask', 'u_mode', 'u_frequency_ref'});
            end
        case 2
            if ~flags.elastic_code
                checkFieldNames(source, {'p0', 'p', 'p_mask', 'p_mode', 'p_frequency_ref', 'ux', 'uy', 'u_mask', 'u_mode', 'u_frequency_ref'});
            else
                checkFieldNames(source, {'p0', 'sxx', 'syy', 'sxy', 's_mask', 's_mode', 'ux', 'uy', 'u_mask', 'u_mode'});
            end
        case 3
            if ~flags.elastic_code
                checkFieldNames(source, {'p0', 'p', 'p_mask', 'p_mode', 'p_frequency_ref', 'ux', 'uy', 'uz', 'u_mask', 'u_mode', 'u_frequency_ref'});
            else
                checkFieldNames(source, {'p0', 'sxx', 'syy', 'szz', 'sxy', 'sxz', 'syz', 's_mask', 's_mode', 'ux', 'uy', 'uz', 'u_mask', 'u_mode'});
            end
    end
    
    % check initial pressure input
    if isfield(source, 'p0')
        
        % set flag
        flags.source_p0 = true;
        
        % check size and contents
        if isempty(source.p0) || ~sum(source.p0(:) ~= 0)
            
            % if the initial pressure is empty, remove field
            source = rmfield(source, 'p0');
            flags.source_p0 = false;
            
        elseif ~all(size(source.p0) == size(kgrid.k))
            
            % throw an error if p0 is not the correct size
            error('source.p0 must be the same size as the computational grid.');
            
        end
        
        % if using the elastic code, reformulate source.p0 in terms of the
        % stress source terms using the fact that source.p = [0.5 0.5] /
        % (2*CFL) is the same as source.p0 = 1 
        if flags.elastic_code
            
            % set flag
            flags.source_p0_elastic = true;
            
            % assign source
            source.s_mask = ones(size(kgrid.k));
            source.sxx(:, 1) = -reshape(source.p0, 1, []) / 2;
            source.sxx(:, 2) = source.sxx(:, 1);
            source.syy(:, 1:2) = repmat(source.sxx(:, 1), [1, 2]);
            if kgrid.dim == 3
                source.szz(:, 1:2) = repmat(source.sxx(:, 1), [1, 2]);
            end
            
        end
    end

    % check for a time varying pressure source input
    if isfield(source, 'p')

        % force p_mask to be given if p is given
        enforceFields(source, {'p_mask'});
        
        % check mask is the correct size
        if (numDim(source.p_mask) ~= kgrid.dim) || (all(size(source.p_mask) ~= size(kgrid.k)))
            error('source.p_mask must be the same size as the computational grid.');
        end
        
        % check mask is not empty
        if sum(source.p_mask(:)) == 0
            error('source.p_mask must be a binary grid with at least one element set to 1.');
        end
        
        % don't allow both source.p0 and source.p in the same simulation
        % USERS: please contact us via http://www.k-wave.org/forum if this
        % is a problem 
        if flags.source_p0
            error('source.p0 and source.p can''t be defined in the same simulation.');
        end
        
        % check the source mode input is valid
        if isfield(source, 'p_mode')
            if ~ischar(source.p_mode) || (~strcmp(source.p_mode, 'additive') && ~strcmp(source.p_mode, 'dirichlet') && ~strcmp(source.p_mode, 'additive-no-correction'))
                error('source.p_mode must be set to ''additive'', ''additive-no-correction'', or ''dirichlet''.');
            end
        else
            source.p_mode = SOURCE_P_MODE_DEF;
        end         
        
        % check if a reference frequency is defined
        if isfield(source, 'p_frequency_ref')
           
            % set flag
            flags.use_w_source_correction_p = true;
            
            % check frequency is a scalar, positive number
            validateattributes(source.p_frequency_ref, {'numeric'}, {'scalar', 'nonnegative'}, '', 'source.p_frequency_ref');
            
            % check frequency is within range
            if source.p_frequency_ref > (kgrid.k_max * min(medium.sound_speed(:)) / (2*pi))
                error('source.p_frequency_ref is higher than the maximum frequency supported by the spatial grid.');
            end
            
            % change source mode to no include k-space correction
            source.p_mode = 'additive-no-correction';
            
        end
        
        % set source flag to the length of the source, this allows source.p
        % to be shorter than kgrid.Nt
        flags.source_p = length(source.p(1, :));
        if flags.source_p > kgrid.Nt
           disp('  WARNING: source.p has more time points than kgrid.Nt, remaining time points will not be used.');
        end        

        % create an indexing variable corresponding to the location of all
        % the source elements 
        p_source_pos_index = find(source.p_mask ~= 0);
        
        % check if the mask is binary or labelled
        p_unique = unique(source.p_mask);
        
        % create a second indexing variable
        if numel(p_unique) <= 2 && sum(p_unique) == 1
            
            % binary source mask
            flags.source_p_labelled = false;
            
            % if more than one time series is given, check the number of time
            % series given matches the number of source elements, or the number
            % of labelled sources
            if (size(source.p, 1) > 1) && (length(source.p(:,1)) ~= sum(source.p_mask(:)))
                error('The number of time series in source.p must match the number of source elements in source.p_mask.');
            end            
            
            % set signal index to all elements
            p_source_sig_index = ':';
            
        else
            
            % labelled source mask
            flags.source_p_labelled = true;
            
            % check the source labels are monotonic, and start from 1
            if (sum(p_unique(2:end) - p_unique(1:end-1)) ~= (numel(p_unique) - 1)) || (~any(p_unique == 1))
                error('If using a labelled source.p_mask, the source labels must be monotonically increasing and start from 1.');
            end
            
            % make sure the correct number of input signals are given
            if size(source.p, 1) ~= (numel(p_unique) - 1)
                error('The number of time series in source.p must match the number of labelled source elements in source.p_mask.');
            end   
            
            % set signal index to the labels (this allows one input signal
            % to be used for each source label)
            p_source_sig_index = source.p_mask(source.p_mask ~= 0);
            
        end
        
        % convert the data type depending on the number of indices
        eval(['p_source_pos_index = ' index_data_type '(p_source_pos_index);']);
        if flags.source_p_labelled
            eval(['p_source_sig_index = ' index_data_type '(p_source_sig_index);']);
        end
                    
    end

    % check for time varying velocity source input and set source flag
    if isfield(source, 'ux') || isfield(source, 'uy') || isfield(source, 'uz') || isfield(source, 'u_mask') 

        % force u_mask to be given
        enforceFields(source, {'u_mask'});
        
        % check mask is the correct size
        if (numDim(source.u_mask) ~= kgrid.dim) || (all(size(source.u_mask) ~= size(kgrid.k)))
            error('source.u_mask must be the same size as the computational grid.');
        end        
        
        % check mask is not empty
        if sum(source.u_mask(:)) == 0
            error('source.u_mask must be a binary grid with at least one element set to 1.');
        end        
        
        % check the source mode input is valid
        if isfield(source, 'u_mode')
            if ~ischar(source.u_mode) || (~strcmp(source.u_mode, 'additive') && ~strcmp(source.u_mode, 'dirichlet') && ~strcmp(source.u_mode, 'additive-no-correction'))
                error('source.u_mode must be set to ''additive'', ''additive-no-correction'', or ''dirichlet''.');
            end
        else
            source.u_mode = SOURCE_U_MODE_DEF;
        end
        
        % check if a reference frequency is defined
        if isfield(source, 'u_frequency_ref')
           
            % set flag
            flags.use_w_source_correction_u = true;
            
            % check frequency is a scalar, positive number
            validateattributes(source.u_frequency_ref, {'numeric'}, {'scalar', 'nonnegative'}, '', 'source.p_frequency_ref');
            
            % check frequency is within range
            if source.u_frequency_ref > (kgrid.k_max * min(medium.sound_speed(:)) / (2*pi))
                error('source.u_frequency_ref is higher than the maximum frequency supported by the spatial grid.');
            end
            
            % change source mode to no include k-space correction
            source.u_mode = 'additive-no-correction';
            
        end
        
        % set source flags to the length of the sources, this allows the
        % inputs to be defined independently and be of any length
        if isfield(source, 'ux')
            flags.source_ux = length(source.ux(1, :));
            if flags.source_ux > kgrid.Nt
               disp('  WARNING: source.ux has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'uy')
            flags.source_uy = length(source.uy(1, :));
            if flags.source_uy > kgrid.Nt
                disp('  WARNING: source.uy has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'uz')
            flags.source_uz = length(source.uz(1, :)) ;
            if flags.source_uz > kgrid.Nt
                disp('  WARNING: source.uz has more time points than kgrid.Nt, remaining time points will not be used.');
            end          
        end    

        % create an indexing variable corresponding to the location of all
        % the source elements 
        u_source_pos_index = find(source.u_mask ~= 0);      
                
        % check if the mask is binary or labelled
        u_unique = unique(source.u_mask);
        
        % create a second indexing variable
        if numel(u_unique) <= 2 && sum(u_unique) == 1
            
            % binary source mask
            flags.source_u_labelled = false;
            
            % if more than one time series is given, check the number of time
            % series given matches the number of source elements
            if (flags.source_ux && (length(source.ux(:,1)) > 1)) || (flags.source_uy && (length(source.uy(:,1)) > 1)) || (flags.source_uz && (length(source.uz(:,1)) > 1))
                if (flags.source_ux && (length(source.ux(:,1)) ~= sum(source.u_mask(:)))) || (flags.source_uy && (length(source.uy(:,1)) ~= sum(source.u_mask(:)))) || (flags.source_uz && (length(source.uz(:,1)) ~= sum(source.u_mask(:))))
                    error('The number of time series in source.ux (etc) must match the number of source elements in source.u_mask.');
                end
            end           
            
            % set signal index to all elements
            u_source_sig_index = ':';
            
        else
            
            % labelled source mask
            flags.source_u_labelled = true;
            
            % check the source labels are monotonic, and start from 1
            if (sum(u_unique(2:end) - u_unique(1:end-1)) ~= (numel(u_unique) - 1)) || (~any(u_unique == 1))
                error('If using a labelled source.u_mask, the source labels must be monotonically increasing and start from 1.');
            end
            
            % if more than one time series is given, check the number of time
            % series given matches the number of source elements
            if (flags.source_ux && (size(source.ux, 1) ~= (numel(u_unique) - 1))) || (flags.source_uy && (size(source.uy, 1) ~= (numel(u_unique) - 1))) || (flags.source_uz && (size(source.uz, 1) ~= (numel(u_unique) - 1)))
                error('The number of time series in source.ux (etc) must match the number of labelled source elements in source.u_mask.');
            end              
            
            % set signal index to the labels (this allows one input signal
            % to be used for each source label)
            u_source_sig_index = source.u_mask(source.u_mask ~= 0);
            
        end        
        
        % convert the data type depending on the number of indices
        eval(['u_source_pos_index = ' index_data_type '(u_source_pos_index);']);   
        if flags.source_u_labelled
            eval(['u_source_sig_index = ' index_data_type '(u_source_sig_index);']);
        end

    end
    
    % check for time varying stress source input and set source flag
    if isfield(source, 'sxx') || isfield(source, 'syy') || isfield(source, 'szz') || ...
       isfield(source, 'sxy') || isfield(source, 'sxz') || isfield(source, 'syz') || ...
       isfield(source, 's_mask') 

        % force s_mask to be given
        enforceFields(source, {'s_mask'});
        
        % check mask is the correct size
        if (numDim(source.s_mask) ~= kgrid.dim) || (all(size(source.s_mask) ~= size(kgrid.k)))
            error('source.s_mask must be the same size as the computational grid.');
        end        
        
        % check mask is not empty
        if sum(source.s_mask(:)) == 0
            error('source.s_mask must be a binary grid with at least one element set to 1.');
        end        
        
        % check the source mode input is valid
        if isfield(source, 's_mode')
            if ~ischar(source.s_mode) || (~strcmp(source.s_mode, 'additive') && ~strcmp(source.s_mode, 'dirichlet'))
                error('source.s_mode must be set to ''additive'' or ''dirichlet''.');
            end
        else
            source.s_mode = SOURCE_S_MODE_DEF;
        end        
        
        % set source flags to the length of the sources, this allows the
        % inputs to be defined independently and be of any length
        if isfield(source, 'sxx')
            flags.source_sxx = length(source.sxx(1, :));
            if flags.source_sxx > kgrid.Nt
               disp('  WARNING: source.sxx has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'syy')
            flags.source_syy = length(source.syy(1, :));
            if flags.source_syy > kgrid.Nt
               disp('  WARNING: source.syy has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'szz')
            flags.source_szz = length(source.szz(1, :));
            if flags.source_szz > kgrid.Nt
               disp('  WARNING: source.szz has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'sxy')
            flags.source_sxy = length(source.sxy(1, :));
            if flags.source_sxy > kgrid.Nt
               disp('  WARNING: source.sxy has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'sxz')
            flags.source_sxz = length(source.sxz(1, :));
            if flags.source_sxz > kgrid.Nt
               disp('  WARNING: source.sxz has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
        if isfield(source, 'syz')
            flags.source_syz = length(source.syz(1, :));
            if flags.source_syz > kgrid.Nt
               disp('  WARNING: source.syz has more time points than kgrid.Nt, remaining time points will not be used.');
            end
        end
 
        % create an indexing variable corresponding to the location of all
        % the source elements 
        s_source_pos_index = find(source.s_mask ~= 0);      
                
        % check if the mask is binary or labelled
        s_unique = unique(source.s_mask);

        % create a second indexing variable
        if numel(s_unique) <= 2 && sum(s_unique) == 1
            
            % binary source mask
            flags.source_s_labelled = false;
            
            % if more than one time series is given, check the number of time
            % series given matches the number of source elements
            if (flags.source_sxx && (length(source.sxx(:,1)) > 1)) || ...
               (flags.source_syy && (length(source.syy(:,1)) > 1)) || ...
               (flags.source_szz && (length(source.szz(:,1)) > 1)) || ...
               (flags.source_sxy && (length(source.sxy(:,1)) > 1)) || ...
               (flags.source_sxz && (length(source.sxz(:,1)) > 1)) || ...
               (flags.source_syz && (length(source.syz(:,1)) > 1))       
                if (flags.source_sxx && (length(source.sxx(:,1)) ~= sum(source.s_mask(:)))) || ...
                   (flags.source_syy && (length(source.syy(:,1)) ~= sum(source.s_mask(:)))) || ...
                   (flags.source_szz && (length(source.szz(:,1)) ~= sum(source.s_mask(:)))) || ...
                   (flags.source_sxy && (length(source.sxy(:,1)) ~= sum(source.s_mask(:)))) || ...
                   (flags.source_sxz && (length(source.sxz(:,1)) ~= sum(source.s_mask(:)))) || ...
                   (flags.source_syz && (length(source.syz(:,1)) ~= sum(source.s_mask(:))))           
                    error('The number of time series in source.sxx (etc) must match the number of source elements in source.s_mask.');
                end
            end            
            
            % set signal index to all elements
            s_source_sig_index = ':';            
            
        else
            
            % labelled source mask
            flags.source_s_labelled = true;
            
            % check the source labels are monotonic, and start from 1
            if (sum(s_unique(2:end) - s_unique(1:end-1)) ~= (numel(s_unique) - 1)) || (~any(s_unique == 1))
                error('If using a labelled source.s_mask, the source labels must be monotonically increasing and start from 1.');
            end
            
            % if more than one time series is given, check the number of time
            % series given matches the number of source elements
            if (flags.source_sxx && (size(source.sxx, 1) ~= (numel(s_unique) - 1))) ||...
               (flags.source_syy && (size(source.syy, 1) ~= (numel(s_unique) - 1))) ||...
               (flags.source_szz && (size(source.szz, 1) ~= (numel(s_unique) - 1))) ||...
               (flags.source_sxy && (size(source.sxy, 1) ~= (numel(s_unique) - 1))) ||...
               (flags.source_sxz && (size(source.sxz, 1) ~= (numel(s_unique) - 1))) ||...
               (flags.source_syz && (size(source.syz, 1) ~= (numel(s_unique) - 1)))
                error('The number of time series in source.sxx (etc) must match the number of labelled source elements in source.u_mask.');
            end              
            
            % set signal index to the labels (this allows one input signal
            % to be used for each source label)
            s_source_sig_index = source.s_mask(source.s_mask ~= 0);            
            
        end
        
        % convert the data type depending on the number of indices
        eval(['s_source_pos_index = ' index_data_type '(s_source_pos_index);']); 
        if flags.source_s_labelled
            eval(['s_source_sig_index = ' index_data_type '(s_source_sig_index);']); 
        end
        
    end
    
    % clean up unused variables
    clear p_unique u_unique s_unique;
    
else
    
    % ----------------------
    % SOURCE IS A TRANSDUCER
    % ----------------------
    
    % if the sensor is a transducer, check that the simulation is in 3D
    if kgrid.dim ~= 3
        error('Transducer inputs are only compatible with 3D simulations.');
    end    
    
    % get the input signal - this is appended with zeros if required to
    % account for the beamforming delays (this will throw an error if the
    % input signal is not defined)
    transducer_input_signal = source.input_signal;
        
    % get the delay mask that accounts for the beamforming delays and
    % elevation focussing; this is used so that a single time series can be
    % applied to the complete transducer mask with different delays
    delay_mask = source.delay_mask;

    % set source flag - this should be the length of signal minus the
    % maximum delay
    flags.transducer_source = length(transducer_input_signal) - max(delay_mask(:));    
    
    % get the active elements mask
    active_elements_mask = source.active_elements_mask;
    
    % get the apodization mask if not set to 'Rectangular' and convert to a
    % linear array
    if ischar(source.transmit_apodization) && strcmp(source.transmit_apodization, 'Rectangular')
        transducer_transmit_apodization = 1;
    else
        transducer_transmit_apodization = source.transmit_apodization_mask;  
        transducer_transmit_apodization = transducer_transmit_apodization(active_elements_mask ~= 0);
    end    
        
    % create indexing variable corresponding to the active elements
    u_source_pos_index = find(active_elements_mask ~= 0);
    
    % convert the data type depending on the number of indices
    eval(['u_source_pos_index = ' index_data_type '(u_source_pos_index);']);     
    
    % convert the delay mask to an indexing variable (this doesn't need to
    % be modified if the grid is expanded) which tells each point in the
    % source mask which point in the input_signal should be used
    delay_mask = delay_mask(active_elements_mask ~= 0);
    
    % convert the data type depending on the maximum value of the delay
    % mask and the length of the source
    max_delay = max(delay_mask(:)) + length(transducer_input_signal) + 1;
    if max_delay < intmax('uint8')
        delay_mask = uint8(delay_mask);
    elseif max_delay < intmax('uint16')
        delay_mask = uint16(delay_mask);
    elseif max_delay < intmax('uint32')
        delay_mask = uint32(delay_mask);               
    end      
    
    % move forward by 1 as a delay of 0 corresponds to the first point in
    % the input signal
    delay_mask = delay_mask + 1;
            
    % clean up unused variables
    clear active_elements_mask;
    
end

% =========================================================================
% CHECK KGRID TIME INPUTS
% =========================================================================

% check kgrid for t_array existance, and create if not defined
if strcmp(kgrid.t_array, 'auto')
    
    % check for time reversal mode
    if flags.time_rev
        error('kgrid.t_array (Nt and dt) must be defined explicitly in time reversal mode.');
    end
    
    % check for time varying sources
    if (~flags.source_p0_elastic) && (...
            flags.source_p   || ...
            flags.source_ux  || flags.source_uy  || flags.source_uz  || ...
            flags.source_sxx || flags.source_syy || flags.source_szz || ...
            flags.source_sxy || flags.source_sxz || flags.source_syz)
        error('kgrid.t_array (Nt and dt) must be defined explicitly when using a time-varying source.');
    end
        
    % create time array
    if flags.elastic_code

        % consider both the shear and compressional sound speeds, but
        % not including shear speeds of zero
        ss = [medium.sound_speed_compression(:); medium.sound_speed_shear(:)];
        ss(ss == 0) = [];

        % create the time array
        if flags.kspace_elastic_code
            kgrid.makeTime(ss, KSPACE_CFL);
        else
            kgrid.makeTime(ss, PSTD_CFL);
        end

        % cleanup unused variables
        clear ss;

    else

        % create the time array using the compressional sound speed
        kgrid.makeTime(medium.sound_speed, KSPACE_CFL);

    end

end
    
% check kgrid.t_array for stability given medium properties
if ~flags.elastic_code

    % calculate the largest timestep for which the model is stable
    dt_stability_limit = checkStability(kgrid, medium);

    % give a warning if the timestep is larger than stability limit allows
    if kgrid.dt > dt_stability_limit
        disp('  WARNING: time step may be too large for a stable simulation.');
    end
    
end

% check for nonuniform grid and assign to flag
flags.nonuniform_grid = kgrid.nonuniform;

% =========================================================================
% CREATE ANONYMOUS FUNCTIONS FOR ALLOCATING VARIABLES
% =========================================================================

% set storage variable type based on data_cast - this enables the
% output variables to be directly created in the data_cast format,
% rather than creating them in double precision and then casting them
switch data_cast
    case 'off'
        castZeros = @(sz) zeros(sz);
        precision = 'double';
    case 'single'
        castZeros = @(sz) zeros(sz, 'single');
        precision = 'single';
    case 'gsingle'
        castZeros = @(sz) gzeros(sz, 'single');
        precision = 'single';
    case 'gdouble'
        castZeros = @(sz) gzeros(sz, 'double');
        precision = 'double';
    case 'gpuArray'
        if strcmp(data_cast_prepend, 'single')
            if verLessThan('distcomp', '6.2')
                
                % original syntax to initialise array of zeros
                castZeros = @(sz) parallel.gpu.GPUArray.zeros(sz, 'single');
                
            else
                
                % current syntax to initialise array of zeros
                castZeros = @(sz) gpuArray.zeros(sz, 'single');
                
            end
            precision = 'single';
        else
            if verLessThan('distcomp', '6.2')
                
                % original syntax to initialise array of zeros
                castZeros = @(sz) parallel.gpu.GPUArray.zeros(sz, 'double');
                
            else
                
                % current syntax to initialise array of zeros
                castZeros = @(sz) gpuArray.zeros(sz, 'double');
                
            end
            precision = 'double';
        end
    case 'kWaveGPUsingle'
        castZeros = @(sz) zeros(sz, GPUsingle);
        precision = 'single';
    case 'kWaveGPUdouble'
        castZeros = @(sz) zeros(sz, GPUdouble);
        precision = 'double';
    otherwise
        error('Unknown ''DataCast'' option');
end

% =========================================================================
% CHECK FOR VALID INPUT COMBINATIONS
% =========================================================================

% enforce density input if velocity sources or output are being used
if ~user_medium_density_input && (flags.source_ux || flags.source_uy || flags.source_uz || flags.record_u || flags.record_u_max || flags.record_u_rms)
    error('medium.density must be explicitly defined if velocity inputs or outputs are used, even in homogeneous media.');
end

% enforce density input if nonlinear equations are being used
if ~user_medium_density_input && flags.nonlinear
    error('medium.density must be explicitly defined if medium.BonA is specified.');
end

% check sensor compatability options for flags.compute_directivity
if flags.use_sensor && kgrid.dim == 2 && flags.compute_directivity && ~flags.binary_sensor_mask && strcmp(cartesian_interp, 'linear')
    error('sensor directivity fields are only compatible with binary sensor masks or ''CartInterp'' set to ''nearest''.');
end

% check for split velocity output
if flags.record_u_split_field && ~flags.binary_sensor_mask
    error('The option sensor.record = {''u_split_field''} is only compatible with a binary sensor mask.');
end

% check input options for data streaming *****
if flags.stream_to_disk
    if (~flags.use_sensor || flags.time_rev)
        error('The optional input ''StreamToDisk'' is currently only compatible with forward simulations using a non-zero sensor mask.');
    elseif isfield(sensor, 'record') && any(ismember(record.flags(2:end), sensor.record)) 
        error('The optional input ''StreamToDisk'' is currently only compatible with sensor.record = {''p''} (the default).');
    end
end

% make sure the PML size is smaller than the grid if PMLInside is true
if flags.pml_inside && ( ...
   (kgrid.dim == 1 && ((pml_x_size*2 > kgrid.Nx))) || ...
   (kgrid.dim == 2 && ~flags.axisymmetric && ((pml_x_size*2 > kgrid.Nx) || (pml_y_size*2 > kgrid.Ny))) || ...
   (kgrid.dim == 2 &&  flags.axisymmetric && ((pml_x_size*2 > kgrid.Nx) || (pml_y_size   > kgrid.Ny))) || ...
   (kgrid.dim == 3 && ((pml_x_size*2 > kgrid.Nx) || (pml_y_size*2 > kgrid.Ny) || (pml_z_size*2 > kgrid.Nz) )) ...
   )
    error('The size of the PML must be smaller than the size of the grid.');
end             

% make sure the PML is inside if using a non-uniform grid
if flags.nonuniform_grid && ~flags.pml_inside
    error('''PMLInside'' must be true for simulations using non-uniform grids.');
end

% check for compatible input options if saving to disk
if (kgrid.dim == 3 && ischar(flags.save_to_disk) && (~flags.use_sensor || ~flags.binary_sensor_mask || flags.time_rev)) 
    error('The optional input ''SaveToDisk'' is currently only compatible with forward simulations using a non-zero binary sensor mask.');
end

% check the record start time is within range
if flags.use_sensor && ( (sensor.record_start_index > kgrid.Nt) || (sensor.record_start_index < 1) )
    error('sensor.record_start_index must be between 1 and the number of time steps.');
end

% ensure 'WSWA' symmetry if using axisymmetric code with 'SaveToDisk'
if flags.axisymmetric && ~strcmp(radial_symmetry, 'WSWA') && ischar(flags.save_to_disk)
    
    % display a warning only if using WSWS symmetry (not WSWA-FFT)
    if strncmp(radial_symmetry, 'WSWS', 4)
        disp('  WARNING: Optional input ''RadialSymmetry'' changed to ''WSWA'' for compatability with ''SaveToDisk''.');
    end
    
    % update setting
    radial_symmetry = 'WSWA';
    
end

% switch off layout plot in time reversal mode
flags.plot_layout = flags.plot_layout && ~flags.time_rev;

% check for automatic plot scaling
if strcmp(plot_scale, 'auto') 
    flags.plot_scale_auto = true;
end

% check for log plot scaling and store the plot scales
if flags.plot_scale_log && ~flags.plot_scale_auto
    alt_plot_scale_lin = plot_scale;
    alt_plot_scale_log = log10(abs(plot_scale) + log_scale_comp_factor) - log10(log_scale_comp_factor);
    alt_plot_scale_log(1) = -alt_plot_scale_log(1);
end

% force visualisation if flags.record_movie is true
if flags.record_movie
    flags.plot_sim = true;
end

% ensure p0 smoothing is switched off if p0 is empty
if ~flags.source_p0
    flags.smooth_p0 = false;
end

% ensure default display mask is switched off if sensor input is empty
if (~flags.use_sensor || flags.blank_sensor) && strcmp(display_mask, 'default')
    display_mask = 'off';
end

% switch off default display mask if using mesh plot
if kgrid.dim == 2 && flags.mesh_plot && strcmp(display_mask, 'default')
    display_mask = 'off';
end

% start log if required
if flags.create_log
    diary([LOG_NAME '.txt']);
end

% update command line status
if flags.time_rev
    disp('  time reversal mode');
end

% check plot scaling if p0 is given
if flags.source_p0 && ~flags.time_rev && flags.plot_sim && ~flags.plot_scale_auto
    
    % find the maximum input pressure amplitude
    if isfield(source, 'p')
        max_val = max([source.p0(:); source.p(:)]);
    else
        max_val = max(source.p0(:));
    end
    
    % check the plot scaling
    if max_val > PLOT_SCALE_WARNING*plot_scale(2) || PLOT_SCALE_WARNING*max_val < plot_scale(2)
        disp('  WARNING: visualisation plot scale may not be optimal for given source.');
    end
    
    clear max_val;
    
end

% cleanup unused variables
clear *_DEF NUM_REQ_INPUT_VARIABLES user_medium_density_input;

% =========================================================================
% UPDATE COMMAND LINE STATUS WITH SIMULATION PARAMETERS
% =========================================================================

% run subscript to display time step, max supported frequency etc.
kspaceFirstOrder_displaySimulationParams;

% =========================================================================
% SMOOTH AND ENLARGE GRIDS
% =========================================================================

% smooth the initial pressure distribution p0 if required, and then restore
% the maximum magnitude
%   NOTE 1: if p0 has any values at the edge of the domain, the smoothing
%   may cause part of p0 to wrap to the other side of the domain
%   NOTE 2: p0 is smoothed before the grid is expanded to ensure that p0 is
%   exactly zero within the PML
%   NOTE 3: for the axisymmetric code, p0 is smoothed assuming WS origin
%   symmetry
if flags.source_p0 && flags.smooth_p0
        
    % update command line status
    disp('  smoothing p0 distribution...');  
    
    if flags.axisymmetric
        switch radial_symmetry
            case {'WSWA-FFT', 'WSWA'}

                % create a new kWave grid object with expanded radial grid
                kgrid_exp = kWaveGrid(kgrid.Nx, kgrid.dx, kgrid.Ny * 4, kgrid.dy);                 
                
                % mirror p0 in radial dimension using WSWA symmetry 
                p0_exp = zeros(kgrid_exp.Nx, kgrid_exp.Ny);
                p0_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1) =         source.p0;
                p0_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2) = -fliplr(source.p0(:, 2:end));
                p0_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3) =        -source.p0;
                p0_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4) =  fliplr(source.p0(:, 2:end));

            case {'WSWS-FFT', 'WSWS'}

                % create a new kWave grid object with expanded radial grid
                kgrid_exp = kWaveGrid(kgrid.Nx, kgrid.dx, kgrid.Ny * 2 - 2, kgrid.dy);                  
                
                % mirror p0 in radial dimension using WSWS symmetry 
                p0_exp = zeros(kgrid_exp.Nx, kgrid_exp.Ny);
                p0_exp(:,            1:kgrid.Ny)       =        source.p0;
                p0_exp(:, kgrid.Ny + 1:kgrid.Ny*2 - 2) = fliplr(source.p0(:, 2:end - 1));
                
        end
        
        % smooth p0
        p0_exp = smooth(p0_exp, true);
        
        % trim back to original size
        source.p0 = p0_exp(:, 1:kgrid.Ny);
        
        % clean up unused variables
        clear kgrid_exp p0_exp;
        
    else
        source.p0 = smooth(source.p0, true);
    end

    % if using the elastic code, reassign the smoothed p0 to the stress
    if flags.elastic_code
        source.sxx(:, 1)  = -reshape(source.p0, 1, []) / 2;
        source.sxx(:, 2)  = source.sxx(:, 1);
        source.syy(:, 1:2) = repmat(source.sxx(:, 1), [1, 2]);
        if kgrid.dim == 3
            source.szz(:, 1:2) = repmat(source.sxx(:, 1), [1, 2]);
        end
    end
    
end
    
% expand the computational grid if the PML is set to be outside the input
% grid defined by the user
if ~flags.pml_inside
    kspaceFirstOrder_expandGridMatrices;
end

% get maximum prime factors
if flags.axisymmetric
    prime_facs = kgrid.highest_prime_factors(radial_symmetry(1:4));
else
    prime_facs = kgrid.highest_prime_factors;
end

% give warning for bad dimension sizes
if max(prime_facs) > HIGHEST_PRIME_FACTOR_WARNING
    prime_facs(prime_facs == 0) = [];
    disp(['  WARNING: Highest prime factors in each dimension are ' num2str(prime_facs)]);
    disp('           Use dimension sizes with lower prime factors to improve speed');
end
clear prime_facs;

% smooth the sound speed distribution if required
if ~flags.elastic_code
    if flags.smooth_c0 && numDim(medium.sound_speed) == kgrid.dim && numel(medium.sound_speed) > 1
        disp('  smoothing sound speed distribution...');  
        medium.sound_speed = smooth(medium.sound_speed);
    end
else
    if flags.smooth_c0
        if numel(medium.sound_speed_compression) > 1
            disp('  smoothing compressional sound speed distribution...');  
            medium.sound_speed_compression = smooth(medium.sound_speed_compression);
        end
        if numel(medium.sound_speed_shear) > 1
            disp('  smoothing shear sound speed distribution...');  
            medium.sound_speed_shear = smooth(medium.sound_speed_shear);     
        end
    end
end
    
% smooth the ambient density distribution if required
if flags.smooth_rho0 && numDim(medium.density) == kgrid.dim && numel(medium.density) > 1
    disp('  smoothing density distribution...');
    medium.density = smooth(medium.density);
end

% =========================================================================
% CREATE SENSOR AND ABSORPTION VARIABLES
% =========================================================================

% define the output variables and mask indices if using the sensor
if flags.use_sensor
    if ~flags.blank_sensor || ischar(flags.save_to_disk)
        if flags.cuboid_corners
            
            % create empty list of sensor indices
            sensor_mask_index = [];
            
            % loop through the list of cuboid corners, and extract the
            % sensor mask indices for each cube                
            for cuboid_index = 1:size(record.cuboid_corners_list, 2)
                
                % create empty binary mask
                temp_mask = false(size(kgrid.k));
                
                % assign cuboid corners to binary mask
                switch kgrid.dim
                    case 1
                        temp_mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(2, cuboid_index)) = 1;
                    case 2
                        temp_mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(3, cuboid_index),...
                                  record.cuboid_corners_list(2, cuboid_index):record.cuboid_corners_list(4, cuboid_index)) = 1;                       
                    case 3
                        temp_mask(record.cuboid_corners_list(1, cuboid_index):record.cuboid_corners_list(4, cuboid_index),...
                                  record.cuboid_corners_list(2, cuboid_index):record.cuboid_corners_list(5, cuboid_index),...
                                  record.cuboid_corners_list(3, cuboid_index):record.cuboid_corners_list(6, cuboid_index)) = 1;
                end
                
                % extract mask indices
                sensor_mask_index = [sensor_mask_index; find(temp_mask ~= 0)]; %#ok<AGROW>
            
            end
            
            % cleanup unused variables
            clear temp_mask
            
        else
        
            % create mask indices (this works for both normal sensor and
            % transducer inputs)
            sensor_mask_index = find(sensor.mask ~= 0);
            
        end

        % convert the data type depending on the number of indices (this saves
        % memory)
        eval(['sensor_mask_index = ' index_data_type '(sensor_mask_index);']); 

    else

        % set the sensor mask index variable to be empty
        sensor_mask_index = [];

    end
end

% run subscript to create storage variables if not saving to disk
if flags.use_sensor && ~ischar(flags.save_to_disk)
    kspaceFirstOrder_createStorageVariables;    
end

% run subscript to create absorption variables for the fluid code based on
% the expanded and smoothed values of the medium parameters (if not saving
% to disk)
if ~flags.elastic_code && ~ischar(flags.save_to_disk)
    kspaceFirstOrder_createAbsorptionVariables;
end

% =========================================================================
% ASSIGN PSEUDONYMS
% =========================================================================

% now that grids have been enlarged and smoothed, shorten commonly used
% field names (these act only as pointers provided that the values aren't
% modified) 
dt = kgrid.dt;
rho0 = medium.density;
if flags.elastic_code
    c0 = medium.sound_speed_compression;
else
    c0 = medium.sound_speed;
end

% =========================================================================
% SCALE SOURCE TERMS
% =========================================================================

% run subscript to scale the source terms based on the expanded and
% smoothed values of the medium parameters
if flags.scale_source_terms
    kspaceFirstOrder_scaleSourceTerms;
end

% =========================================================================
% CREATE PML INDICES
% =========================================================================
        
% define index variables to remove the PML from the display if the optional
% input 'PlotPML' is set to false
if ~flags.plot_pml
    switch kgrid.dim
        case 1
            x1 = pml_x_size + 1; 
            x2 = kgrid.Nx - pml_x_size;
        case 2
            x1 = pml_x_size + 1;
            x2 = kgrid.Nx - pml_x_size;
            if flags.axisymmetric
                y1 = 1;
            else
                y1 = pml_y_size + 1;
            end
            y2 = kgrid.Ny - pml_y_size;                
        case 3
            x1 = pml_x_size + 1;
            x2 = kgrid.Nx - pml_x_size;
            y1 = pml_y_size + 1;
            y2 = kgrid.Ny - pml_y_size;    
            z1 = pml_z_size + 1;
            z2 = kgrid.Nz - pml_z_size;            
    end
else
    switch kgrid.dim
        case 1
            x1 = 1;
            x2 = kgrid.Nx;
        case 2
            x1 = 1;
            x2 = kgrid.Nx;
            y1 = 1;
            y2 = kgrid.Ny;
        case 3
            x1 = 1;
            x2 = kgrid.Nx;
            y1 = 1;
            y2 = kgrid.Ny;
            z1 = 1;
            z2 = kgrid.Nz;            
    end
end    

% define index variables to allow original grid size to be maintained for
% the _final and _all output variables if 'PMLInside' is set to false
if ~flags.pml_inside
    switch kgrid.dim
        case 1
            record.x1_inside = pml_x_size + 1;
            record.x2_inside = kgrid.Nx - pml_x_size; 
        case 2
            record.x1_inside = pml_x_size + 1;
            record.x2_inside = kgrid.Nx - pml_x_size;
            if flags.axisymmetric
                record.y1_inside = 1;
            else
                record.y1_inside = pml_y_size + 1;
            end
            record.y2_inside = kgrid.Ny - pml_y_size;    
        case 3
            record.x1_inside = pml_x_size + 1;
            record.x2_inside = kgrid.Nx - pml_x_size;
            record.y1_inside = pml_y_size + 1;
            record.y2_inside = kgrid.Ny - pml_y_size;    
            record.z1_inside = pml_z_size + 1;
            record.z2_inside = kgrid.Nz - pml_z_size;    
    end
else
    switch kgrid.dim
        case 1
            record.x1_inside = 1;
            record.x2_inside = kgrid.Nx;
        case 2
            record.x1_inside = 1;
            record.x2_inside = kgrid.Nx;
            record.y1_inside = 1;
            record.y2_inside = kgrid.Ny;
        case 3
            record.x1_inside = 1;
            record.x2_inside = kgrid.Nx;
            record.y1_inside = 1;
            record.y2_inside = kgrid.Ny;
            record.z1_inside = 1;
            record.z2_inside = kgrid.Nz;            
    end
end