%KWAVETRANSDUCER Class definition for k-Wave linear array transducer.
%
% DESCRIPTION:
%     kWaveTransducer creates an object of the kWaveTransducer class which
%     can be substituted for the source or sensor inputs when using 
%     kspaceFirstOrder3D.  
%
%     Note: This function will not work with older versions of MATLAB in
%     which custom class definitions are not supported. 
%
% USAGE:
%     transducer = makeTransducer(kgrid, settings)
%
% INPUTS:
%       
%     kgrid           - k-Wave grid object returned by kWaveGrid
%     settings        - input structure used to define the properties of
%                       the transducer (see below) 
%
%     The following parameters can be appended as fields to the settings
%     input structure. These parameters are fixed when the transducer is
%     initialised (they cannot be changed without re-creating the
%     transducer). All the settings are optional and are given their
%     default values if not defined.    
%
%     number_elements - total number of transducer elements (default = 128)
%     element_width   - width of each element in grid points (default = 1)
%     element_length  - length of each element in grid points
%                       (default = 10)
%     element_spacing - spacing (kerf width) between the transducer
%                       elements in grid points (default = 0)
%     position        - position of the corner of the transducer within the
%                       grid in grid points (default = [1, 1, 1])
%     radius          - radius of curvature of the transducer [m]
%                       (currently only inf is supported) (default = inf) 
%     input_signal    - signal used to drive the ultrasound transducer
%                       (where the sampling rate is defined by kgrid.dt
%                       or kgrid.t_array) (default = [])
%
%     The following parameters can also be appended as fields to the
%     settings input structure, and can additionally be modified after
%     the transducer has been initialised.
% 
%     active_elements - transducer elements that are currently active
%                       elements (default = all elements)
%     beamforming_delay_offset 
%                     - beamforming delay offset (used to force beamforming
%                       delays to be positive) (default = 'auto')  
%     elevation_focus_distance 
%                     - focus depth in the elevation direction [m]
%                       (default = inf)  
%     focus_distance  - focus distance used to calculate beamforming
%                       delays [m] (default = inf)  
%     receive_apodization 
%                     - receive apodization; can be set to any of the
%                       window shapes supported by getWin, or given as a
%                       vector the same length as number_elements
%                       (default = 'Rectangular')   
%     sound_speed     - sound speed used to calculate beamforming delays
%                       [m/s] (default = 1540) 
%     steering_angle  - steering angle used to calculate beamforming delays
%                       [deg] (default = 0)
%     transmit_apodization 
%                     - transmit apodization; can be set to any of the
%                       window shapes supported by getWin, or given as a
%                       vector the same length as number_elements 
%                       (default = 'Rectangular')  
%
% OUTPUTS:
%
%     transducer      - kWaveTransducer object which can be used to replace
%                       the source or sensor inputs of kspaceFirstOrder3D  
% 
%     In addition to the input parameters given above (which are also
%     accessible after the transducer has been created) the kWaveTransducer
%     object has a number of dependent properties and methods.   
%
%     input_signal    - user defined input signal appended and prepended
%                       with additional zeros depending on the values of
%                       focus_distance, elevation_focus_distance, and
%                       steering_angle 
%     mask            - binary mask of the active transducer elements
%     number_active_elements 
%                     - current number of active transducer elements  
%     transducer_width 
%                     - total width of the transducer in grid points
%
%     active_elements_mask 
%                     - return a binary mask of the active transducer
%                       elements (identical to mask)  
%     all_elements_mask 
%                     - return a binary mask of all the transducer
%                       elements (both active and inactive) 
%     beamforming_delays 
%                     - return a vector of the beam forming delays (in
%                       units of time samples) for each active element
%                       based on the focus and steering angle settings 
%     combine_sensor_data(sensor_data)
%                     - combine the sensor data returned by
%                       kspaceFirstOrder-OMP and kspaceFirstOrder-CUDA
%                       to give a single time series per active transducer
%                       element (anlogous to the output from the MATLAB
%                       code), rather than a time series per grid point 
%     delay_mask      - return a mask of the active transducer elements,
%                       where the mask values contain the beamforming
%                       delays (an integer input can also be given to
%                       control the beamforming delays used, where 1:
%                       both delays, 2: elevation only, 3: azimuth only)
%     elevation_beamforming_delays 
%                     - return a vector of the elevation beam forming
%                       delays (in units of time samples) for each active
%                       element based on the elevation focus setting
%     get_receive_apodization 
%                     - return the receive apodization
%     get_transmit_apodization 
%                     - return the transmit apodization
%     indexed_active_elements_mask 
%                     - return a mask of the active transducer elements,
%                       where the mask values indicate which transducer
%                       element each grid point corresponds to   
%     indexed_elements_mask 
%                     - return a mask of all the transducer elements (both
%                       active and inactive), where the mask values 
%                       indicate which transducer element each grid point
%                       corresponds to  
%     plot            - plot the transducer using voxelPlot
%     properties      - print a list of the transducer properties to the
%                       command line 
%     scan_line(sensor_data)
%                     - combine the input sensor data using the current
%                       apodization and beamforming setting to generate a
%                       single scan line  
%     transmit_apodization_mask 
%                     - return a mask of the active transducer elements,
%                       where the mask values contain the apodization
%                       weights
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 9th December 2010
%     last update     - 11th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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

% Define as a handle class - A handle class constructor returns a handle
% object that is a reference to the object created. You can assign the
% handle object to multiple variables or pass it to functions without
% causing MATLAB to make a copy of the original object. A function that
% modifies a handle object passed as an input argument does not need to
% return the object.
classdef kWaveTransducer < handle

    % define the properties of the transducer that cannot be modified by
    % the user after the transducer is initialised (these parameters are
    % stored). The numbers assigned here are the default values used if the
    % parameters are not set explicitly by the user.
    properties (GetAccess = 'public', SetAccess = 'private')
        
        % the total number of transducer elements
        number_elements = 128;        
        
        % the width of each element in grid points
        element_width = 1;
        
        % the length of each element in grid points
        element_length = 20;
        
        % the spacing (kerf width) between the transducer elements in grid
        % points
        element_spacing = 0;    
        
        % the position of the corner of the transducer in the grid
        position = [1, 1, 1];      
        
        % the radius of curvature of the transducer [m]
        radius = inf;        
        
    end
    
    % define the properties of the transducer that can be modified by the
    % user after the transducer is initialised (these parameters are
    % stored). The numbers assigned here are the default values used if the
    % parameters are not set explicitly by the user.
    properties (Access = 'public')
        
        % the transducer elements that are currently active elements
        active_elements;        
               
        % the focus depth in the elevation direction [m]
        elevation_focus_distance = inf;

        % transmit apodization
        transmit_apodization = 'Rectangular';        
        
        % receive apodization
        receive_apodization = 'Rectangular';
        
        % sound speed used to calculate beamforming delays [m/s]
        sound_speed = 1540;
        
        % focus distance used to calculate beamforming delays [m]
        focus_distance = inf;
        
        % time index to start recording if transducer is used as a sensor
        record_start_index = 1;
        
    end
    
    % define the dependent properties (these parameters are computed when
    % queried). Some of these properties act like normal user-modifiable
    % properties (e.g., steering_angle), but instead have explicit set and
    % get functions that change an internal stored value. This allows error
    % checking when the properties are changed.
    properties(Dependent = true)
        
        % binary mask of the active transducer elements
        active_elements_mask;
        
        % indexed mask of the active transducer elements, where the index
        % indicates the transducer element that each grid point belongs to.
        indexed_active_elements_mask;
        
        % binary mask of all the transducer elements (both active and
        % inactive)
        all_elements_mask;
        
        % indexed mask of all the transducer elements, where the index
        % indicates the transducer element that each grid point belongs to.
        indexed_elements_mask;
        
        % pseudonym for the active_elements_mask to maintain compatability
        % when the transducer is used in place of sensor or source.
        mask;
        
        % total width of the transducer in grid points
        transducer_width;
        
        % current number of active transducer elements
        number_active_elements;
        
        % user defined input signal with appended zeros dependent on the
        % beamforming (steering_angle, focus_distance) settings
        input_signal;
        
        % maximum steering angle - this sets the beamforming delay offset
        steering_angle_max;
        
        % number of zeros appended to input signal to allow a single time
        % series to be used within kspaceFirstOrder3D (either set to 'auto'
        % or based on the setting for steering_angle_max)
        appended_zeros;
        
        % offset used to make all the delays in the delay_mask positive
        % (either set to 'auto' or based on the setting for
        % steering_angle_max) 
        beamforming_delays_offset;
        
        % steering angle used to calculate beamforming delays [deg]
        steering_angle;
        
    end
    
    % define the hidden properties - these cannot be seen or accessed
    % directly by the user, and are used to store internal properties
    properties (Hidden = true, Access = 'private')
        
        % indexed sensor mask from which the other masks are computed
        indexed_mask;
        
        % indexed mask for the voxels in each element used to compute
        % elevation focussing
        indexed_element_voxel_mask;
        
        % size of the grid in which the transducer is defined
        stored_grid_size;
        
        % corresponding grid spacing
        grid_spacing; 
        
        % time spacing
        dt;

        % original copy of the user defined input signal
        stored_input_signal;
        
        % stored value of appended_zeros (accessed using get and
        % set methods). This is used to set the number of zeros that are
        % appended and prepended to the input signal.
        stored_appended_zeros = 'auto';
        
        % stored value of the minimum beamforming delay. This is used to
        % offset the delay mask so that all the delays are >= 0
        stored_beamforming_delays_offset = 'auto';
        
        % stored value of the steering_angle_max (accessed using get and
        % set methods). This can be set by the user and is used to derive
        % the two parameters above.
        stored_steering_angle_max = 'auto';
               
        % stored value of the steering_angle (accessed using get and
        % set methods)
        stored_steering_angle = 0;
        
    end
    
    % constructor function
    methods
        function transducer = kWaveTransducer(kgrid, transducer_properties)
            
            % allocate the grid size and spacing
            transducer.stored_grid_size = [kgrid.Nx, kgrid.Ny, kgrid.Nz];
            transducer.grid_spacing = [kgrid.dx, kgrid.dy, kgrid.dz];
            
            % allocate the temporal spacing
            if isnumeric(kgrid.dt)
                transducer.dt = kgrid.dt;
            elseif isnumeric(kgrid.t_array)
                transducer.dt = kgrid.t_array(2) - kgrid.t_array(1);
            else
                error('kgrid.dt or kgrid.t_array must be explicitly defined');
            end
            
            % check for properties input
            if nargin == 1
                transducer_properties = [];
            end
            
            % check the input fields
            checkFieldNames(transducer_properties, {...
                'number_elements', 'element_width', 'element_length', 'element_spacing',...
                'position', 'radius', ...
                'active_elements', 'elevation_focus_distance', 'focus_distance', 'steering_angle',...
                'receive_apodization', 'transmit_apodization',...
                'sound_speed', 'input_signal', 'steering_angle_max'});
            
            % replace default settings with user defined properties if
            % given
            % -------------------------------------------------------------
            
            if isfield(transducer_properties, 'number_elements')
                
                % force value to be a positive integer
                transducer.number_elements = round(abs(transducer_properties.number_elements));
                
            end
    
            if isfield(transducer_properties, 'element_width')
                
                % force value to be a positive integer
                transducer.element_width = round(abs(transducer_properties.element_width));
                
            end
            
            if isfield(transducer_properties, 'element_length')
                
                % force value to be a positive integer
                transducer.element_length = round(abs(transducer_properties.element_length));
                
            end
            
            if isfield(transducer_properties, 'element_spacing')
                
                % force value to be a positive integer
                transducer.element_spacing = round(abs(transducer_properties.element_spacing));
                
            end
            
            if isfield(transducer_properties, 'position')
                
                % force values to be positive integers
                transducer.position = round(abs(transducer_properties.position));
                
            end
            
            if isfield(transducer_properties, 'radius')
                transducer.radius = transducer_properties.radius;
                
                % only allow an infinite radius for now
                if ~isinf(transducer.radius)
                    error('Only a value of transducer.radius = inf is currently supported');
                end
            end
            
            if isfield(transducer_properties, 'active_elements')
                transducer.active_elements = transducer_properties.active_elements;
            else
                transducer.active_elements = ones(transducer_properties.number_elements, 1);
            end
            
            if isfield(transducer_properties, 'elevation_focus_distance')
                transducer.elevation_focus_distance = transducer_properties.elevation_focus_distance;
            end
            
            if isfield(transducer_properties, 'receive_apodization')
                
                % if a user defined apodization is given explicitly, check
                % the length of the input
                if isnumeric(transducer_properties.receive_apodization) ...
                        && (length(transducer_properties.receive_apodization) ~= transducer.number_active_elements)
                    error('The length of the receive apodization input must match the number of active elements');
                end
                                       
                % assign the input                  
                transducer.receive_apodization = transducer_properties.receive_apodization;    
                
            end
            
            if isfield(transducer_properties, 'transmit_apodization')
                
                % if a user defined apodization is given explicitly, check
                % the length of the input
                if isnumeric(transducer_properties.transmit_apodization) ...
                        && (length(transducer_properties.transmit_apodization) ~= transducer.number_active_elements)
                    error('The length of the transmit apodization input must match the number of active elements');
                end
                    
                % assign the input
                transducer.transmit_apodization = transducer_properties.transmit_apodization; 
                
            end            
                                
            if isfield(transducer_properties, 'sound_speed')
                
                transducer.sound_speed = transducer_properties.sound_speed;
                
                % check to see the sound_speed is positive
                if ~(transducer.sound_speed> 0)
                    error('transducer.sound_speed must be greater than 0');
                end    
                
            end
            
            if isfield(transducer_properties, 'focus_distance')
                transducer.focus_distance = transducer_properties.focus_distance;
            end
            
            if isfield(transducer_properties, 'input_signal')
                
                transducer.stored_input_signal = transducer_properties.input_signal;
                
                % force the input signal to be a column vector
                if numDim(transducer.input_signal) == 1
                    transducer.stored_input_signal = reshape(transducer.stored_input_signal, [], 1);
                else
                    error('transducer.input_signal must be a one-dimensional array');
                end
                
            end
            
            if isfield(transducer_properties, 'steering_angle_max')
                
                % set the maximum steering angle using the set method (this
                % avoids having to duplicate error checking)
                transducer.steering_angle_max = transducer_properties.steering_angle_max;
                
            end
            
            if isfield(transducer_properties, 'steering_angle')
                
                % set the maximum steering angle using the set method (this
                % avoids having to duplicate error checking)
                transducer.steering_angle = transducer_properties.steering_angle;
                
            end
            
            % -------------------------------------------------------------
            
            % check the transducer fits into the grid
            if sum(transducer.position(:) == 0)
                error('The defined transducer must be positioned within the grid');
            elseif (transducer.position(2) + transducer.number_elements*transducer.element_width + (transducer.number_elements - 1)*transducer.element_spacing) > transducer.stored_grid_size(2)
                error('The defined transducer is too large or positioned outside the grid in the y-direction');
            elseif (transducer.position(3) + transducer.element_length) > transducer.stored_grid_size(3)
                transducer.position(3)
                transducer.element_length
                transducer.stored_grid_size(3)
                error('The defined transducer is too large or positioned outside the grid in the z-direction');
            elseif transducer.position(1) > transducer.stored_grid_size(1)
                error('The defined transducer is positioned outside the grid in the x-direction');
            end
         
            % assign the data type for the transducer matrix based on the
            % number of different elements (uint8 supports 255 numbers so
            % most of the time this data type will be used)
            if transducer.number_elements < intmax('uint8')
                mask_type = 'uint8';
            elseif transducer.number_elements < intmax('uint16')
                mask_type = 'uint16';
            elseif transducer.number_elements < intmax('uint32')
                mask_type = 'uint32';                
            else
                mask_type = 'uint64';
            end
            
            % create an empty transducer mask (the grid points within
            % element 'n' are all given the value 'n')
            transducer.indexed_mask = zeros(transducer.stored_grid_size, mask_type);
                
            % create a second empty mask used for the elevation beamforming
            % delays (the grid points across each element are numbered 1 to
            % M, where M is the number of grid points in the elevation
            % direction)
            transducer.indexed_element_voxel_mask = zeros(transducer.stored_grid_size, mask_type);
                
            % create the corresponding indexing variable 1:M
            element_voxel_index = repmat(1:transducer.element_length, transducer.element_width, 1);
                            
            % for each transducer element, calculate the grid point indices
            for element_index = 1:transducer.number_elements

                % assign the current transducer position
                element_pos_x = transducer.position(1);
                element_pos_y = transducer.position(2) + (transducer.element_width + transducer.element_spacing)*(element_index - 1);
                element_pos_z = transducer.position(3);
                
                % assign the grid points within the current element the
                % index of the element
                transducer.indexed_mask(element_pos_x, element_pos_y:element_pos_y + transducer.element_width - 1, element_pos_z:element_pos_z + transducer.element_length - 1) = element_index;
                
                % assign the individual grid points an index corresponding
                % to their order across the element
                transducer.indexed_element_voxel_mask(element_pos_x, element_pos_y:element_pos_y + transducer.element_width - 1, element_pos_z:element_pos_z + transducer.element_length - 1) = element_voxel_index;

            end

            % double check the transducer fits within the desired grid size
            if any(size(transducer.indexed_elements_mask) ~= transducer.stored_grid_size)
                error('Desired transducer is larger than the input grid_size');
            end
            
        end
    end
    
    % set and get functions for dependent variables
    methods
        
        % return input signal
        function signal = get.input_signal(obj)
            
            % copy the stored input signal
            signal = obj.stored_input_signal;
            
            % check the signal is not empty
            if isempty(signal)
                error('Transducer input signal is not defined');
            end
            
            % automatically prepend and append zeros if the beamforming
            % delay offset is set
            
            % check if the beamforming delay offset is set. If so, use this
            % number to prepend and append this number of zeros to the
            % input signal. Otherwise, calculate how many zeros are needed
            % and prepend and append these.
            if ~strcmp(obj.stored_appended_zeros, 'auto')
                
                % use the current value of the beamforming offset to add
                % zeros to the input signal
                signal = [zeros(obj.stored_appended_zeros, 1); signal; zeros(obj.stored_appended_zeros, 1)];
                
            else
            
                % get the current delay beam forming 
                delay_mask = obj.delay_mask;

                % find the maximum delay
                delay_max = max(delay_mask(:));

                % count the number of leading zeros in the input signal
                leading_zeros = find(signal ~= 0, 1) - 1;

                % count the number of trailing zeros in the input signal
                trailing_zeros = find(flipud(signal) ~= 0, 1) - 1;

                % check the number of leading zeros is sufficient given the
                % maximum delay
                if leading_zeros < delay_max + 1
                    
                    disp(['  prepending transducer.input_signal with ' num2str(delay_max - leading_zeros + 1) ' leading zeros']);

                    % prepend extra leading zeros
                    signal = [zeros(delay_max - leading_zeros + 1, 1); signal];  
                    
                end

                % check the number of leading zeros is sufficient given the
                % maximum delay
                if trailing_zeros < delay_max + 1
                    
                    disp(['  appending transducer.input_signal with ' num2str(delay_max - trailing_zeros + 1) ' trailing zeros']);

                    % append extra trailing zeros
                    signal = [signal; zeros(delay_max - trailing_zeros + 1, 1)];
                    
                end
            end
            
        end
        
        function transducer_width = get.transducer_width(obj)
            
            % return the overall length of the transducer
            transducer_width = obj.number_elements * obj.element_width + (obj.number_elements - 1) * obj.element_spacing;
            
        end
        
        function number_active_elements = get.number_active_elements(obj)
            
            % return the current number of active elements
            number_active_elements = sum(obj.active_elements(:));
            
        end
        
        % allow mask query to allow compatability with regular sensor
        % structure - return the active sensor mask
        function mask = get.mask(obj)
            mask = obj.active_elements_mask;
        end    
        
        % return a binary mask showing the locations of the active elements
        function mask = get.active_elements_mask(obj)
            
            % copy the indexed elements mask
            mask = obj.indexed_mask;
            
            % remove inactive elements from the mask
            for element_index = 1:obj.number_elements
                mask(mask == element_index) = obj.active_elements(element_index);
            end
            
            % convert remaining mask to binary
            mask(mask ~= 0) = 1;

        end      
        
        % return a binary mask showing the locations of all the elements
        % (both active and inactive)
        function mask = get.all_elements_mask(obj)
            mask = obj.indexed_mask;
            mask(mask ~= 0) = 1;
        end        
        
        % return a copy of the indexed elements mask
        function mask = get.indexed_elements_mask(obj)
            mask = obj.indexed_mask;
        end
        
        % return a copy of the indexed elements mask only for the active
        % elements
        function mask = get.indexed_active_elements_mask(obj)
            
            % copy the indexed elements mask
            mask = obj.indexed_mask;
            
            % remove inactive elements from the mask
            for element_index = 1:obj.number_elements
                if ~(obj.active_elements(element_index))
                    mask(mask == element_index) = 0;
                end
            end
            
            % force the lowest element index to be 1
            lowest_active_element_index = find(obj.active_elements, 1);
            mask(mask ~= 0) = mask(mask ~= 0) - lowest_active_element_index + 1;

        end   
     
        % return the stored value of the steering angle
        function steering_angle = get.steering_angle(obj)
            steering_angle = obj.stored_steering_angle;
        end
        
        % set the stored value of the steering angle
        function set.steering_angle(obj, steering_angle)
            
            % force to be scalar
            steering_angle = steering_angle(1);
            
            % check if the steering angle is between -90 and 90
            if ~(steering_angle < 90 && steering_angle > -90)
                error('Input for steering_angle must be betweeb -90 and 90 degrees.');
            end
            
            % check if the steering angle is less than the maximum steering
            % angle
            if ~strcmp(obj.stored_steering_angle_max, 'auto') && (abs(steering_angle) > obj.stored_steering_angle_max)
                error('Input for steering_angle cannot be greater than steering_angle_max.');
            end
            
            % update the stored value
            obj.stored_steering_angle = steering_angle;
            
        end 
        
        % return the stored value of the maximum steering angle
        function steering_angle_max = get.steering_angle_max(obj)
            steering_angle_max = obj.stored_steering_angle_max;
        end
        
        % set the stored value of the maximum steering angle and
        % beamforming offset
        function set.steering_angle_max(obj, steering_angle_max)
            
            % force input to be scalar and positive
            steering_angle_max = steering_angle_max(1);
            
            % check the steering angle is within range
            if ~(steering_angle_max < 90 && steering_angle_max > -90)
                error('Input for steering_angle_max must be between -90 and 90.');
            end
            
            % check the maximum steering angle is greater than the current
            % steering angle
            if ~strcmp(obj.stored_steering_angle_max, 'auto') && (abs(obj.stored_steering_angle) > steering_angle_max)
                error('Input for steering_angle_max cannot be less than the current steering_angle.');
            end
            
            % overwrite the stored value
            obj.stored_steering_angle_max = steering_angle_max;
            
            % store a copy of the current value for the steering angle
            current_steering_angle = obj.stored_steering_angle;
            
            % overwrite with the user defined maximum value
            obj.stored_steering_angle = steering_angle_max;
            
            % set the beamforming delay offset to zero (this means the
            % delays will contain negative values which we can use to
            % derive the new values for the offset)
            obj.stored_appended_zeros = 0;
            obj.stored_beamforming_delays_offset = 0;
            
            % get the element beamforming delays and reverse
            delay_times = -obj.beamforming_delays;            
            
            % get the maximum and minimum beamforming delays
            max_delay = max(delay_times);
            min_delay = min(delay_times);
            
            % add the maximum and minimum elevation delay if the elevation
            % focus is not set to infinite 
            if ~isinf(obj.elevation_focus_distance)
                max_delay = max_delay + max(obj.elevation_beamforming_delays);
                min_delay = min_delay + min(obj.elevation_beamforming_delays); 
            end
            
            % set the beamforming offset to the difference between the
            % maximum and minimum delay 
            obj.stored_appended_zeros = max_delay - min_delay;
            
            % set the minimum beamforming delay (converting to a positive
            % number)
            obj.stored_beamforming_delays_offset = - min_delay;
            
            % reset the previous value of the steering angle
            obj.stored_steering_angle = current_steering_angle;
            
        end
        
        % return the stored value of appended zeros
        function appended_zeros = get.appended_zeros(obj)
            appended_zeros = obj.stored_appended_zeros;
        end
        
        % return the stored value of the offset used to force the values in
        % delay_mask to be >= 0
        function beamforming_delays_offset = get.beamforming_delays_offset(obj)
            beamforming_delays_offset = obj.stored_beamforming_delays_offset;
        end
        
    end
    
    % general class methods
    methods
        
        % calculate the beamforming delays based on the focus and steering
        % settings
        function delay_times = beamforming_delays(obj)
            
            % calculate the element pitch in [m]
            element_pitch = (obj.element_spacing + obj.element_width) * obj.grid_spacing(2);

            % create indexing variable
            element_index = -(obj.number_active_elements - 1)/2:(obj.number_active_elements - 1)/2;           

            % check for focus depth
            if isinf(obj.focus_distance)
                
                % calculate time delays for a steered beam  
                delay_times = element_pitch * element_index * sin(obj.steering_angle * pi / 180) / (obj.sound_speed);        % [s]
                
            else   
                
                % calculate time delays for a steered and focussed beam
                delay_times = obj.focus_distance / obj.sound_speed * (1 - sqrt(...
                    1 ... 
                    + (element_index * element_pitch ./ obj.focus_distance).^2 ... 
                    - 2 * (element_index * element_pitch ./ obj.focus_distance) * sin(obj.steering_angle * pi / 180) ...
                    ));      % [s]
                
            end
                        
            % convert the delays to be in units of time points
            delay_times = round(delay_times / obj.dt);            
            
        end
        
        % calculate the elevation beamforming delays based on the focus 
        % setting
        function delay_times = elevation_beamforming_delays(obj)
            
            if ~isinf(obj.elevation_focus_distance)
            
                % create indexing variable
                element_index = -(obj.element_length - 1)/2:(obj.element_length - 1)/2; 

                % calculate time delays for a focussed beam
                delay_times = (obj.elevation_focus_distance - sqrt( (element_index * obj.grid_spacing(3)).^2 + obj.elevation_focus_distance^2 )) ./ obj.sound_speed;

                % convert the delays to be in units of time points and then
                % reverse 
                delay_times = -round(delay_times / obj.dt);
                
            else
               
                % create an empty array
                delay_times = zeros(1, obj.element_length);
                
            end

        end
        
        % function to create a scan line based on the input sensor data and
        % the current apodization and beamforming setting
        function line = scan_line(obj, sensor_data)
           
            % get the current apodization setting
            apodization = obj.get_receive_apodization;
            
            % get the current beamforming weights and reverse
            delays = -obj.beamforming_delays;         

            % offset the received sensor_data by the beamforming delays and
            % apply receive apodization 
            for element_index = 1:obj.number_active_elements
                if delays(element_index) > 0
                    
                    % shift element data forwards
                    sensor_data(element_index, :) = apodization(element_index).*[sensor_data(element_index, 1 + delays(element_index):end), zeros(1, delays(element_index))];
                    
                elseif delays(element_index) < 0
                    
                    % shift element data backwards
                    sensor_data(element_index, :) = apodization(element_index).*[zeros(1, -delays(element_index)), sensor_data(element_index, 1:end + delays(element_index))];
                    
                end
            end

            % form the a-line summing across the elements
            line = sum(sensor_data);
                       
        end
        
        % function to average the individual time series recorded at each
        % grid point of each transducer element (as returned by the C++
        % code) and return a single time series per active transducer
        % element (as returned by the MATLAB code)
        function sensor_data_sum = combine_sensor_data(obj, sensor_data)
            
            % check the data is the correct size
            if size(sensor_data, 1) ~= (obj.number_active_elements * obj.element_width * obj.element_length)
                error('The number of time series in the input sensor_data must match the number of grid points in the active tranducer elements.');
            end
            
            % get index of which element each time series belongs to  
            ind = obj.indexed_active_elements_mask(obj.indexed_active_elements_mask > 0);

            % create empty output
            sensor_data_sum = zeros(obj.number_active_elements, size(sensor_data, 2));
            
            % check if an elevation focus is set
            if isinf(obj.elevation_focus_distance)
                
                % loop over time series and sum
                for ii = 1:length(ind)
                    sensor_data_sum(ind(ii), :) = sensor_data_sum(ind(ii), :) + sensor_data(ii, :);
                end
                    
            else
                
                % get the elevation delay for each grid point of each
                % active transducer element (this is given in units of grid
                % points)
                dm = obj.delay_mask(2);
                dm = dm(obj.active_elements_mask ~= 0);
                
                % loop over time series, shift and sum
                for ii = 1:length(ind)
                    sensor_data_sum(ind(ii), 1:end - dm(ii)) = sensor_data_sum(ind(ii), 1:end - dm(ii)) + sensor_data(ii, 1 + dm(ii):end);
                end
                
            end
            
            % divide by number of time series in each element
            sensor_data_sum = sensor_data_sum .* (1 / (obj.element_width * obj.element_length));
            
        end
        
        % plot the transducer using voxelPlot
        function plot(obj)
            voxelPlot(double(obj.all_elements_mask));
        end
        
        % print out the transducer properties
        function properties(obj)
            disp(' ');
            disp('k-Wave Transducer Properties');
            disp(['  transducer position:       [' num2str(obj.position) ']']);
            disp(['  transducer width:          ' scaleSI(obj.transducer_width*obj.grid_spacing(2)) 'm (' num2str(obj.transducer_width) ' grid points)']);
            disp(['  number of elements:        ' num2str(obj.number_elements)]);
            disp(['  number of active elements: ' num2str(obj.number_active_elements) ' (elements ' num2str(find(obj.active_elements, 1, 'first')) ' to ' num2str(find(obj.active_elements, 1, 'last')) ')']);
            disp(['  element width:             ' scaleSI(obj.element_width*obj.grid_spacing(2)) 'm (' num2str(obj.element_width) ' grid points)']);
            disp(['  element spacing (kerf):    ' scaleSI(obj.element_spacing*obj.grid_spacing(2)) 'm (' num2str(obj.element_spacing) ' grid points)']);
            disp(['  element pitch:             ' scaleSI((obj.element_spacing + obj.element_width)*obj.grid_spacing(2)) 'm (' num2str((obj.element_spacing + obj.element_width)) ' grid points)']);
            disp(['  element length:            ' scaleSI(obj.element_length*obj.grid_spacing(3)) 'm (' num2str(obj.element_length) ' grid points)']);
            disp(['  sound speed:               ' num2str(obj.sound_speed) 'm/s']);
            if isinf(obj.focus_distance)
                disp('  focus distance:             infinite');
            else
                disp(['  focus distance:            ' scaleSI(obj.focus_distance) 'm']);
            end
            if isinf(obj.elevation_focus_distance)
                disp('  elevation focus distance:  infinite');
            else
                disp(['  elevation focus distance:  ' scaleSI(obj.elevation_focus_distance) 'm']);
            end
            disp(['  steering angle:            ' num2str(obj.steering_angle) ' degrees']);
            disp(['  steering angle max:        ' num2str(obj.steering_angle_max) ' degrees']);
        end
      
    end
    
    % hidden class methods used with kspaceFirstOrder3D
    methods (Hidden = true, Access = 'public') 
        
        % return a delay mask based on the transmit beamforming delays.
        % This is used if the transducer is a source. The delays are set
        % for each grid point of the active_elements_mask and define which
        % element of the input_signal should be used at each source point.
        % This mask is incremented within kspaceFirstOrder3D after each
        % time step. The mask is indexed from 0 and is incremented by 1
        % within kspaceFirstOrder3D .
        function mask = delay_mask(obj, mode)
            % mode == 1: both delays
            % mode == 2: elevation only
            % mode == 3: azimuth only
                                  
            % assign the delays to a new mask using the indexed_element_mask
            indexed_active_elements_mask_copy = obj.indexed_active_elements_mask;
            active_elements_index = find(indexed_active_elements_mask_copy ~= 0);
            mask = zeros(obj.stored_grid_size);

            % calculate azimuth focus delay times provided they are not all
            % zero
            if (~isinf(obj.focus_distance) || (obj.steering_angle ~= 0)) && (nargin == 1 || mode ~= 2)
                
                % get the element beamforming delays and reverse
                delay_times = -obj.beamforming_delays;
                
                % add delay times
                mask(active_elements_index) = delay_times(indexed_active_elements_mask_copy(active_elements_index));
                
            end

            % calculate elevation focus time delays provided each element
            % is longer than one grid point
            if ~isinf(obj.elevation_focus_distance) && (obj.element_length > 1) && (nargin == 1 || mode ~= 3)
                
                % get elevation beamforming delays
                elevation_delay_times = obj.elevation_beamforming_delays;
                
                % get current mask
                element_voxel_mask = obj.indexed_element_voxel_mask;

                % add delay times
                mask(active_elements_index) = mask(active_elements_index) + elevation_delay_times(element_voxel_mask(active_elements_index)).';
                
            end
            
            % shift delay times (these should all be >= 0, where a value of
            % 0 means no delay)
            if strcmp(obj.stored_appended_zeros, 'auto')
                mask(active_elements_index) = mask(active_elements_index) - min(min(min(mask(active_elements_index)))); 
            else
                mask(active_elements_index) = mask(active_elements_index) + obj.stored_beamforming_delays_offset;
            end            
                        
        end
        
        % return a binary mask based on the elevation focus. This sets the
        % elements to extract at each timestep from the sensor_data_buffer
        % used in kspaceFirstOrder3D. This is used if the transducer is a
        % sensor and the elevation focus is set. The mask is of the form:  
        %
        %   1 0 0 0
        %   0 1 0 0
        %   0 0 1 0
        %   0 0 1 0
        %   0 0 0 1
        %   0 0 0 1
        %   0 0 0 1
        %   0 0 1 0
        %   0 0 1 0
        %   0 1 0 0
        %   1 0 0 0
        %
        % The mask remains constant, but extracts new elements at each time
        % step as the FIFO buffer sensor_data_buffer is updated
        function mask = elevation_beamforming_mask(obj)
            
            % get elevation beamforming mask
            delay_mask = obj.delay_mask(2);
            
            % extract the active elements
            delay_mask = delay_mask(obj.active_elements_mask ~= 0);
            
            % force delays to start from zero
            delay_mask = delay_mask - min(delay_mask(:));
            
            % create an empty output mask
            mask = zeros(length(delay_mask), max(delay_mask) + 1);
            
            % populate the mask by setting 1's at the index given by the
            % delay time
            for index = 1:length(delay_mask)
                mask(index, delay_mask(index) + 1) = 1;
            end
            
            % flip the mask so the shortest delays are at the right
            mask = fliplr(mask);
            
        end
        
        % convert the transmit apodization into the form of a element mask,
        % where the apodization values are placed at the grid points
        % belonging to the active transducer elements. These values are
        % then extracted in the correct order within
        % kspaceFirstOrder_inputChecking using apodization =
        % transmit_apodization_mask(active_elements_mask ~= 0) 
        function mask = transmit_apodization_mask(obj)
            
            % get transmit apodization
            apodization = obj.get_transmit_apodization;
            
            % create an empty mask;
            mask = zeros(obj.stored_grid_size);
            
            % assign the apodization values to every grid point in the
            % transducer
            mask_index = obj.indexed_active_elements_mask;
            mask_index = mask_index(mask_index ~= 0);
            mask(obj.active_elements_mask == 1) = apodization(mask_index);
            
        end
        
        % allow the mask to be resized by the simulation functions (this is
        % used when 'PMLInside' is set to false)
        function expand_grid(obj, expand_size)
            obj.indexed_mask = expandMatrix(obj.indexed_mask, expand_size, 0);
        end
        
        % allow the mask to be resized by the simulation functions (this is
        % used when 'PMLInside' is set to false)
        function retract_grid(obj, retract_size)
            obj.indexed_mask = obj.indexed_mask(1 + retract_size(1):end - retract_size(1), 1 + retract_size(2):end - retract_size(2), 1 + retract_size(3):end - retract_size(3));
        end 
        
        % return the size of the grid
        function grid_size = grid_size(obj)
            grid_size = obj.stored_grid_size;
        end
        
    end
    
    % internal class methods only accessible by other functions within
    % kWaveTransducer 
    methods (Hidden = true, Access = 'private') 
        
        % return the transmit apodization, converting strings of window
        % type to actual numbers using getWin
        function apodization = get_transmit_apodization(obj)
            
            % check if a user defined apodization is given and whether this
            % is still the correct size (in case the number of active
            % elements has changed)
            if isnumeric(obj.transmit_apodization)
                if (length(obj.transmit_apodization) ~= obj.number_active_elements)
                    error('The length of the transmit apodization input must match the number of active elements');
                else
                    
                    % assign apodization
                    apodization = obj.transmit_apodization;
                    
                end
            else
                
                % if the number of active elements is greater than 1,
                % create apodization using getWin, otherwise, assign 1
                if obj.number_active_elements > 1
                    apodization = getWin(obj.number_active_elements, obj.transmit_apodization);
                else
                    apodization = 1;
                end
                
            end            
                     
        end                
        
        % return the receive apodization, converting strings of window
        % type to actual numbers using getWin
        function apodization = get_receive_apodization(obj)
            
            % check if a user defined apodization is given and whether this
            % is still the correct size (in case the number of active
            % elements has changed)
            if isnumeric(obj.receive_apodization)
                if (length(obj.receive_apodization) ~= obj.number_active_elements)
                    error('The length of the receive apodization input must match the number of active elements');
                else
                    
                    % assign apodization
                    apodization = obj.receive_apodization;
                    
                end
            else
                
                % if the number of active elements is greater than 1,
                % create apodization using getWin, otherwise, assign 1
                if obj.number_active_elements > 1
                    apodization = getWin(obj.number_active_elements, obj.receive_apodization);
                else
                    apodization = 1;
                end
                
            end  
            
        end    
        
    end
end