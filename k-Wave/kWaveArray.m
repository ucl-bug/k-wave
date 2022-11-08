%KWAVEARRAY Class definition for k-Wave array.
%
% DESCRIPTION:
%     kWaveArray is a class definition that allows the definition of
%     multi-element transducer arrays in 1D, 2D, and 3D for use with the
%     k-Wave simulation functions. The individual elements are defined
%     using physical (rather than grid) parameters in Cartesian
%     coordinates, and their representation on the grid is automatically
%     calculated using off-grid sources as described in [1]. 
%
%     There are two key advantages of using the kWaveArray class to define
%     arrays. (1) The geometry is stair-case free. (2) The transmit signals
%     (for sources) and receive signals (for sensors) are defined per
%     physical transducer element, rather than per grid point.
%
%     No information about the actual simulation grid is stored within
%     objects of the kWaveArray class. The idea is that an array can be
%     defined using physical characteristics, and then re-used for many
%     different simulations. Note, some methods do use the grid information
%     for calculations, and therefore take an object of the kWaveGrid class
%     as an input. The implementation assumes that the grid sampling is
%     uniform, i.e., dx = dy = dz.
%
%     To calculate the off-grid source, first a uniform sampling of
%     Cartesian points is generated covering the element shape (e.g., bowl,
%     disc, etc). For each sample point, a band-limited interpolant is
%     computed corresponding to a point source at that location. These
%     point source responses are summed and scaled to give the source mask
%     for that element. This means the source mask will be non-local, i.e.,
%     will extend spatially beyond where the physical source lies on the
%     grid. Note, any Cartesian points that are outside the grid are
%     automatically removed. 
%
%     In the current version, objects of the kWaveArray class cannot be
%     passed directly to the simulation functions in place of the source
%     and sensor structures. Instead, methods of the kWaveArray class can
%     be used to automatically create these inputs and process the outputs.
%     See examples for further details.
%
%     [1] Wise, E. S., Cox, B. T., Jaros, J., & Treeby, B. E. (2019).
%     Representing arbitrary acoustic source and sensor distributions in
%     Fourier collocation methods. The Journal of the Acoustical Society of
%     America, 146(1), 278-288.
%
% USAGE:
%     karray = kWaveArray
%     karray = kWaveArray(...)
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the
%     default computational settings. See [1] for more details.
%
%     'BLITolerance'    - Scalar value controlling where the spatial extent
%                         of the BLI at each point is trunctated as a
%                         portion of the maximum value (default = 0.05).
%     'BLIType'         - String controlling the BLI expression that is
%                         used for each point source, either 'sinc' or
%                         'exact' (default = 'sinc'). BLITolerance is
%                         ignored if 'exact' is specified.
%     'SinglePrecision' - Boolean controlling whether the calculation of
%                         the off-grid mask and distributed source are
%                         performed in single precision to improve
%                         performance.  
%     'UpsamplingRate'  - Oversampling used to distribute the off-grid
%                         points compared to the equivalent number of
%                         on-grid points (default = 10).
%
% OUTPUTS:
%     karray            - kWaveArray object which can be used to define
%                         source and sensor inputs for the kWave
%                         simulation functions.
%
% PROPERTIES:
%     Properties of the kWaveArray class which can be queried, but not
%     directly modified.
%
%     .array_transformation 
%                       - Position of the array defined by an affine
%                         transform. Used to move the array without needing
%                         to redefine the individual positions of the
%                         elements. Set by using the .setArrayPosition and
%                         .setAffineTransform methods.
%     .dim              - Number of spatial dimensions.
%     .elements         - Structure containing the properties for each
%                         element in the array. New elements are added
%                         using the .addXXXElement methods.
%     .number_elements  - Number of transducer elements in the array.
%
% METHODS:
%     Methods of the kWaveArray class. After creating a kWaveArray object,
%     methods are called using the syntax:
%
%         karray.methodName(inputs)
%         output = karray.methodName(inputs)
%
%     The first syntax is used for methods with no output, and the second
%     for methods with an output.
% 
%     addAnnularArray(position, radius, diameters, focus_pos)
%
%         Add annular array (3D simulations). The array elements are
%         indexed from the inner-most element outwards. If radius ~= inf,
%         then the annular array is positioned over the surface of a bowl.
%
%         position  - Centre of rear surface [bx, by, bz] [m].
%         radius    - Radius of curvature of a bowl on which the annular
%                     array lies [m]. 
%         diameters - 2 x num_elements array containing pairs of inner and
%                     outer aperture diameter (diameter of opening) [m]. 
%         focus_pos - Any point on the beam axis [fx, fy, fz] [m].
%
%     addAnnularElement(position, radius, diameters, focus_pos)
%
%         Add annular element (3D simulations). Setting diameter = [0,
%         diameter] will generate a bowl the same as addBowlElement.
%
%         position  - Centre of rear surface [bx, by, bz] [m].
%         radius    - Radius of curvature of a bowl on which the annular
%                     element lies [m]. 
%         diameters - 2 element array containing inner and outer aperture
%                     diameter (diameter of opening) [m].  
%         focus_pos - Any point on the beam axis [fx, fy, fz] [m].
%
%     addArcElement(position, radius, diameter, focus_pos)
%
%         Add arc-shaped element to the array (2D simulations).
%
%         position  - Centre of rear surface (arc midpoint) [bx, by] [m].
%         radius    - Radius of curvature of arc [m].
%         diameter  - Diameter of arc opening [m].
%         focus_pos - Any point on the beam axis [fx, fy] [m].
%
%     addBowlElement(position, radius, diameter, focus_pos)
%
%         Add bowl-shaped element to the array (3D simulations).
%
%         position  - Centre of rear surface [bx, by, bz] [m].
%         radius    - Radius of curvature of bowl [m].
%         diameter  - Diameter of bowl opening [m].
%         focus_pos - Any point on the beam axis [fx, fy, fz] [m].
%     
%     addCustomElement(integration_points, measure, element_dim, label)
%
%         Add custom volume, area, or line element to the array defined as
%         a series of integration points (see [1]; 1D/2D/3D simulations).
%         Note, for custom elements, the integration point density should
%         be sufficient relative the intended grid spacing (see [1]).
%
%         integration_points - 1 x num_points (1D), 2 x num_points (2D) or
%                              3 x num_points (3D) array of Cartesian
%                              coordinates [m].
%         measure            - Length [m] of line element, area [m^2] of
%                              area element, or volume [m^3] of volume
%                              element corresponding to element_dim.
%         element_dim        - Integer specifying whether the custom
%                              element is a line element (1), area element
%                              (2), or volume element (3). Note,
%                              element_dim does not need to match the
%                              dimension of the simulation, i.e., it is
%                              possible to use a line element in a 3D
%                              array.
%         label              - String identifier.
%
%     addDiscElement(position, diameter, focus_pos)
%
%         Add disc-shaped element to the array (2D/3D simulations).
%
%         position  - Centre of disc surface [bx, by] or [bx, by, bz] [m].
%         diameter  - Diameter of the disc [m].
%         focus_pos - Any point on beam axis [fx, fy, fz] (not used for 2D
%                     simulations) [m]. 
%
%     addRectElement(position, Lx, Ly, theta)
%
%         Add rectangular element to the array (2D/3D simulations). The
%         rectangle is created in the x-y plane and then rotated.
%
%         position - Centre of rect [bx, by] or [bx, by, bz] [m].
%         Lx       - Height of rect (along x-axis before rotation) [m].
%         Ly       - Width of rect (along y-axis before rotation) [m].
%         theta    - Either a scalar (2D) or three element vector (3D) 
%                    [tx, ty, tz] specifying the orientation of the
%                    rectangle [deg]. In 3D, the rotations are specified
%                    about x-y'-z'' (intrinsic rotations) or z-y-x
%                    (extrinsic rotations). All rotations are
%                    counter-clockwise. Can be set to [] if no rotation.
%
%     addLineElement(start_point, end_point)
%
%         Add line element to the array (1D/2D/3D simulations).
%
%         start_point - Start coordinate for the line given as a one (1D),
%                       two (2D), or three (3D) element vector [m]. 
%         end_point   - End coordinate for the line given as a one (1D),
%                       two (2D), or three (3D) element vector [m]. 
%
%     combined_sensor_data = combineSensorData(kgrid, sensor_data)
%
%         When using karray.getArrayBinaryMask to define a sensor mask for
%         the k-Wave simulation functions, the returned sensor data is
%         defined for every grid point that forms part of the array. This
%         method combines the sensor data with the appropriate weights and
%         returns a single time series (or value) for each physical array
%         element (rather than each grid point). The data is returned in
%         the same order as the transducer elements were added to the
%         array. 
%
%         combined_sensor_data - Combined sensor data (one time series per
%                                transducer element).
%         kgrid                - Grid object returned by kWaveGrid.
%         sensor_data          - Sensor data returned by simulation
%                                functions for a sensor mask given by
%                                getArrayBinaryMask. 
%
%     mask = getArrayBinaryMask(kgrid)
%
%         Returns binary mask containing all grid points that form part of
%         the array. Note, as the array elements use off-grid sources (see
%         [1]), this will be non-local to some extent (depending on the
%         values set for 'BLITolerance' and 'BLIType').
%
%         mask  - Binary mask (matrix of 1s and 0s) specifying the grid
%                 points that form part of the array.
%         kgrid - Grid object returned by kWaveGrid.
%
%     grid_weights = getArrayGridWeights(kgrid)
%   
%         Returns matrix containing sum of off-grid source weights for each
%         grid point in the domain (defined by kgrid).
%
%         grid_weights - Matrix of grid weights.
%         kgrid        - Grid object returned by kWaveGrid.
%
%     distributed_source = getDistributedSourceSignal(kgrid, source_signal)
%
%         When defining a source input, a binary source mask (e.g,
%         source.p_mask) should be defined using getArrayBinaryMask, and
%         the time varying source (e.g., source.p) should be defined using
%         this method. This automatically calculates the appropriate
%         weighted source signal for each grid point that forms part of the
%         off-grid source. The source data is assigned in the same order as
%         the transducer elements were added to the array.
%
%         distributed_source - Source signal used to define time-varying
%                              source inputs for the k-Wave simulations
%                              functions, e.g, source.p. 
%         kgrid              - Grid object returned by kWaveGrid.
%         source_signal      - Source signal for each transducer element
%                              defined as an array [number_elements, Nt].
%
%     mask = getElementBinaryMask(kgrid, element_num)
%
%         Returns binary mask containing all grid points that form part of
%         the specified element. Note, as the array elements use off-grid
%         sources (see [1]), this will be non-local (depending on the
%         values set for 'BLITolerance' and 'BLIType').
%
%         mask         - Binary mask (matrix of 1s and 0s) specifying the
%                        grid points that form part of the element.
%         kgrid        - Grid object returned by kWaveGrid.
%         element_num  - Element number to return binary mask for.
%
%     grid_weights = getElementGridWeights(kgrid, element_num)
%
%         Returns matrix containing off-grid source weights for the
%         transducer element defined as element_num (see [1]).
%
%         grid_weights - Matrix of grid weights.
%         kgrid        - Grid object returned by kWaveGrid.
%         element_num  - Element number to return grid weights for.
%
%     element_pos = getElementPositions
%
%         Returns a [dim, num_elements] array containing the positions of
%         each array element after applying array_transformation (if
%         defined).
%
%         element_pos - Matrix of element positions [m].
%
%     plotArray(new_figure)
%
%         Plot the array elements. If new_figure is true (the default), a
%         new figure window is created, otherwise the elements are added to
%         the currently active figure window.
%
%         new_figure - Boolean controlling whether a new figure window is
%                      created.
%
%     removeElement(element_num)
%
%         Remove specified element from the array.
%
%         element_num  - Element number to remove from the array.
%
%     setAffineTransform(affine_transform)
%
%         Sets value for array_transformation defined as an affine
%         transformation matrix.
%
%         affine_transform - Affine transform given as a [3, 3] matrix in
%                            2D or [4, 4] matrix in 3D.
%
%     setArrayPosition(translation, rotation)
%
%         Sets the property array_transformation (an affine transform)
%         based on the values for translation and rotation. The
%         translations are given as [dx, dy] in 2D and [dx, dy, dz] in 3D.
%         The rotations angle/s are given as [th] in 2D (counter-clockwise)
%         and [x_th, y_th, z_th] in 3D (rotation about x then y' then z'').
% 
%         translation - Array translation [m].
%         rotation    - Array rotation [degrees].
%
%     setOptionalInputs(...)
%
%         Method to define the optional inputs (see OPTIONAL INPUTS above).
%         Can be called when the kWaveArray object is defined, or later.
%
% ABOUT:
%     author      - Bradley Treeby and Elliott Wise
%     date        - 5th September 2018
%     last update - 31st July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2021 Bradley Treeby and Elliott Wise
%
% See also offGridPoints, getDeltaBLI

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
classdef kWaveArray < handle

    % define the private properties. The numbers assigned here are the default values used if
    % the parameters are not set by the constructor.
    properties (GetAccess = public, SetAccess = private)
        
        % structure holding the properties of each element
        elements                = {};
        
        % integer number of elements in the array
        number_elements         = 0;
        
        % number of group_id (used to link sub-groups of elements when
        % adding arrays instead of elements)
        number_groups           = 0;        
        
        % number of dimensions (2 for 2D, 3 for 3D)
        dim                     = 0;
        
        % position and orientation of the whole array defined by an affine
        % transformation matrix
        array_transformation    = [];
        
        % flag for axisymmetric coordinate system
        axisymmetric            = false;
        
    end
    
    % define internal properties (these can be accessed but not altered by
    % the user after the object is created, the parameters are also hidden
    % from view) 
    properties (Hidden = true, GetAccess = public, SetAccess = private)
        
        % off-grid parameters
        bli_tolerance    = 0.05;
        bli_type         = 'sinc';
        upsampling_rate  = 10;
        single_precision = false;
        
    end     
    
    % define literals
    properties (Constant, Hidden = true)
    
        % option to use spiral sampling for disc elements
        use_spiral_disc_points = false;
        
        % number of points used to plot arcs
        num_arc_plot_points = 100;
        
        % element plot colour (BUG blue)
        element_plot_colour = [0, 158, 194]/255;
        
    end
    
    % constructor function
    methods
        function karray = kWaveArray(varargin)
            
            % replace defaults with user defined values if provided
            if nargin > 0
                karray.setOptionalInputs(varargin{:});
            end
            
        end
    end
    
    % general methods
    methods
        
        % ----------------------------------
        % ADD AND REMOVE TRANSDUCER ELEMENTS
        % ----------------------------------
        
        % add annular array
        function addAnnularArray(obj, position, radius, diameters, focus_pos)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real', 'numel', 3},       'addAnnularArray', 'position',  1);
            validateattributes(radius,    {'numeric'}, {'finite', 'real', 'numel', 1},       'addAnnularArray', 'radius',    2);
            validateattributes(diameters, {'numeric'}, {'finite', 'real', 'size', [2, NaN]}, 'addAnnularArray', 'diameters', 3);
            validateattributes(focus_pos, {'numeric'}, {'finite', 'real', 'numel', 3},       'addAnnularArray', 'focus_pos', 4);
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = 3;
            end
            
            % check that the element is being added to an array in 3D
            if obj.dim ~= 3
                error(['3D annular array cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end
            
            % number of elements in the annular array
            annular_array_num_el = size(diameters, 2);
            
            % increment the number of groups
            obj.number_groups = obj.number_groups + 1;
            
            % add the individual elements
            for el_ind = 1:annular_array_num_el

                % increment the number of elements
                obj.number_elements = obj.number_elements + 1;
                
                % compute arc angle from chord: https://en.wikipedia.org/wiki/Chord_(geometry)
                varphi_min = asin(diameters(1, el_ind) ./ (2 * radius));
                varphi_max = asin(diameters(2, el_ind) ./ (2 * radius));

                % compute area
                area = 2 * pi .* radius.^2 .* (1 - cos(varphi_max)) - 2 * pi .* radius.^2 .* (1 - cos(varphi_min));

                % add the element properties
                obj.elements{obj.number_elements}.group_id            = obj.number_groups;
                obj.elements{obj.number_elements}.group_type          = 'annular_array';
                obj.elements{obj.number_elements}.element_number      = el_ind;
                obj.elements{obj.number_elements}.type                = 'annulus';
                obj.elements{obj.number_elements}.dim                 = 2;
                obj.elements{obj.number_elements}.position            = position(:).';
                obj.elements{obj.number_elements}.radius_of_curvature = radius;
                obj.elements{obj.number_elements}.inner_diameter      = diameters(1, el_ind);
                obj.elements{obj.number_elements}.outer_diameter      = diameters(2, el_ind);
                obj.elements{obj.number_elements}.focus_position      = focus_pos(:).';
                obj.elements{obj.number_elements}.active              = true;
                obj.elements{obj.number_elements}.measure             = area;
                
            end

        end
        
        % add annular array
        function addAnnularElement(obj, position, radius, diameters, focus_pos)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real', 'numel', 3}, 'addAnnularElement', 'position',  1);
            validateattributes(radius,    {'numeric'}, {'finite', 'real', 'numel', 1}, 'addAnnularElement', 'radius',    2);
            validateattributes(diameters, {'numeric'}, {'finite', 'real', 'numel', 2}, 'addAnnularElement', 'diameters', 3);
            validateattributes(focus_pos, {'numeric'}, {'finite', 'real', 'numel', 3}, 'addAnnularElement', 'focus_pos', 4);
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = 3;
            end
            
            % check that the element is being added to an array in 3D
            if obj.dim ~= 3
                error(['3D annular array cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;

            % compute arc angle from chord: https://en.wikipedia.org/wiki/Chord_(geometry)
            varphi_min = asin(diameters(1) ./ (2 * radius));
            varphi_max = asin(diameters(2) ./ (2 * radius));

            % compute area
            area = 2 * pi .* radius.^2 .* (1 - cos(varphi_max)) - 2 * pi .* radius.^2 .* (1 - cos(varphi_min));

            % add the element properties
            obj.elements{obj.number_elements}.group_id            = 0;
            obj.elements{obj.number_elements}.type                = 'annulus';
            obj.elements{obj.number_elements}.dim                 = 2;
            obj.elements{obj.number_elements}.position            = position(:).';
            obj.elements{obj.number_elements}.radius_of_curvature = radius;
            obj.elements{obj.number_elements}.inner_diameter      = diameters(1);
            obj.elements{obj.number_elements}.outer_diameter      = diameters(2);
            obj.elements{obj.number_elements}.focus_position      = focus_pos(:).';
            obj.elements{obj.number_elements}.active              = true;
            obj.elements{obj.number_elements}.measure             = area;

        end
        
        % add arc element
        function addArcElement(obj, position, radius, diameter, focus_pos)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real', 'numel', 2}, 'addArcElement', 'position',  1);
            validateattributes(radius,    {'numeric'}, {'finite', 'real', 'numel', 1}, 'addArcElement', 'radius',    2);
            validateattributes(diameter,  {'numeric'}, {'finite', 'real', 'numel', 1}, 'addArcElement', 'diameter',  3);
            validateattributes(focus_pos, {'numeric'}, {'finite', 'real', 'numel', 2}, 'addArcElement', 'focus_pos', 4);
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = 2;
            end
            
            % check that the element is being added to an array in 2D
            if obj.dim ~= 2
                error(['2D arc element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end            
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % compute arc angle from chord: https://en.wikipedia.org/wiki/Chord_(geometry)
            varphi_max = asin(diameter ./ (2 * radius));
            
            % compute length
            length = 2 * radius * varphi_max;
            
            % add the element properties
            obj.elements{obj.number_elements}.group_id            = 0;
            obj.elements{obj.number_elements}.type                = 'arc';
            obj.elements{obj.number_elements}.dim                 = 1;
            obj.elements{obj.number_elements}.position            = position(:).';
            obj.elements{obj.number_elements}.radius_of_curvature = radius;
            obj.elements{obj.number_elements}.diameter            = diameter;
            obj.elements{obj.number_elements}.focus_position      = focus_pos(:).';
            obj.elements{obj.number_elements}.active              = true;
            obj.elements{obj.number_elements}.measure             = length;

        end
        
        % add bowl element
        function addBowlElement(obj, position, radius, diameter, focus_pos)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real', 'numel', 3}, 'addBowlElement', 'position',  1);
            validateattributes(radius,    {'numeric'}, {'finite', 'real', 'numel', 1}, 'addBowlElement', 'radius',    2);
            validateattributes(diameter,  {'numeric'}, {'finite', 'real', 'numel', 1}, 'addBowlElement', 'diameter',  3);
            validateattributes(focus_pos, {'numeric'}, {'finite', 'real', 'numel', 3}, 'addBowlElement', 'focus_pos', 4);
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = 3;
            end
            
            % check that the element is being added to an array in 3D
            if obj.dim ~= 3
                error(['3D bowl element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % compute arc angle from chord: https://en.wikipedia.org/wiki/Chord_(geometry)
            varphi_max = asin(diameter ./ (2 * radius));
            
            % compute area
            area = 2 * pi .* radius.^2 .* (1 - cos(varphi_max));
            
            % add the element properties
            obj.elements{obj.number_elements}.group_id            = 0;
            obj.elements{obj.number_elements}.type                = 'bowl';
            obj.elements{obj.number_elements}.dim                 = 2;
            obj.elements{obj.number_elements}.position            = position(:).';
            obj.elements{obj.number_elements}.radius_of_curvature = radius;
            obj.elements{obj.number_elements}.diameter            = diameter;
            obj.elements{obj.number_elements}.focus_position      = focus_pos(:).';
            obj.elements{obj.number_elements}.active              = true;
            obj.elements{obj.number_elements}.measure             = area;

        end
        
        % add custom element, where integration_points is a 1 x num_points,
        % 2 x num_points or 3 x num_points array of Cartesian coordinates
        function addCustomElement(obj, integration_points, measure, element_dim, label)
            
            % check inputs
            validateattributes(integration_points, {'numeric'}, {'finite', 'real', '2d'},                             'addCustomElement', 'integration_points', 1); 
            validateattributes(measure,            {'numeric'}, {'finite', 'real', 'numel', 1},                       'addCustomElement', 'measure',            2);
            validateattributes(element_dim,        {'numeric'}, {'finite', 'integer', 'numel', 1,  '>=', 1, '<=', 3}, 'addCustomElement', 'element_dim',        3);
            validateattributes(label,              {'char'},    {},                                                   'addCustomElement', 'label',              4);
            
            % check the dimensionality of the integration points
            input_dim = size(integration_points, 1);
            if (input_dim < 1) || (input_dim > 3)
                error('Input integration_points must be a 1 x N (in 1D), 2 x N (in 2D), or 3 x N (in 3D) array.');
            end
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = input_dim;
            end
            
            % check that the element is being added to an array with the
            % correct dimensions
            if obj.dim ~= input_dim
                error([num2str(input_dim) 'D custom element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end            
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % store the integration points
            obj.elements{obj.number_elements}.group_id           = 0;
            obj.elements{obj.number_elements}.type               = 'custom';
            obj.elements{obj.number_elements}.dim                = element_dim;
            obj.elements{obj.number_elements}.label              = label;
            obj.elements{obj.number_elements}.integration_points = integration_points;
            obj.elements{obj.number_elements}.active             = true;
            obj.elements{obj.number_elements}.measure            = measure;
            
        end
        
        % add disc element
        function addDiscElement(obj, position, diameter, focus_pos)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real'},             'addDiscElement', 'position',  1);
            validateattributes(diameter,  {'numeric'}, {'finite', 'real', 'numel', 1}, 'addDiscElement', 'diameter',  2);
            
            % check if element is in 2D or 3D
            if numel(position) == 2
                coord_dim = 2;
            elseif numel(position) == 3
                coord_dim = 3;
            else
                error('Input position for disc element must be specified as a 2 (2D) or 3 (3D) element array.');
            end
            
            % check focus input if in 3D
            if coord_dim == 3
                validateattributes(focus_pos, {'numeric'}, {'finite', 'real', 'numel', 3}, 'addDiscElement', 'focus_pos', 3);
            else
                focus_pos = [];
            end
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = coord_dim;
            end
            
            % check element is being added to array of correct dimensions
            if obj.dim ~= coord_dim
                error([num2str(coord_dim) 'D disc element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % compute area
            area = pi .* (diameter / 2).^2;
            
            % add the element properties
            obj.elements{obj.number_elements}.group_id       = 0;
            obj.elements{obj.number_elements}.type           = 'disc';
            obj.elements{obj.number_elements}.dim            = 2;
            obj.elements{obj.number_elements}.position       = position(:).';
            obj.elements{obj.number_elements}.diameter       = diameter;
            obj.elements{obj.number_elements}.focus_position = focus_pos(:).';
            obj.elements{obj.number_elements}.active         = true;
            obj.elements{obj.number_elements}.measure        = area;

        end
        
        % add rectangular element
        function addRectElement(obj, position, Lx, Ly, theta)
            
            % check inputs
            validateattributes(position,  {'numeric'}, {'finite', 'real'},             'addRectElement', 'position',  1);
            validateattributes(Lx,        {'numeric'}, {'finite', 'real', 'numel', 1}, 'addRectElement', 'Lx',        2);
            validateattributes(Ly,        {'numeric'}, {'finite', 'real', 'numel', 1}, 'addRectElement', 'Ly',        3);
            
            % check if element is in 2D or 3D
            if numel(position) == 2
                coord_dim = 2;
                theta_length = 1;
            elseif numel(position) == 3
                coord_dim = 3;
                theta_length = 3;
            else
                error('Input position for rectangular element must be specified as a 2 (2D) or 3 (3D) element array.');
            end
            
            % check angle input
            validateattributes(theta, {'numeric'}, {'finite', 'real', 'numel', theta_length}, 'addRectElement', 'theta', 4);
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = coord_dim;
            end
            
            % check element is being added to array of correct dimensions
            if obj.dim ~= coord_dim
                error([num2str(coord_dim) 'D rectangular element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % compute area
            area = Lx .* Ly;
            
            % add the element properties
            obj.elements{obj.number_elements}.group_id    = 0;
            obj.elements{obj.number_elements}.type        = 'rect';
            obj.elements{obj.number_elements}.dim         = 2;
            obj.elements{obj.number_elements}.position    = position(:).';
            obj.elements{obj.number_elements}.length      = Lx;
            obj.elements{obj.number_elements}.width       = Ly;
            obj.elements{obj.number_elements}.orientation = theta(:).';
            obj.elements{obj.number_elements}.active      = true;
            obj.elements{obj.number_elements}.measure     = area;
            
        end
        
        % add line element
        function addLineElement(obj, start_point, end_point)
            
            % check inputs
            validateattributes(start_point, {'numeric'}, {'finite', 'real'}, 'addLineElement', 'start_point', 1);
            validateattributes(end_point,   {'numeric'}, {'finite', 'real'}, 'addLineElement', 'end_point',   2);
            
            % check the dimensionality of the start and end points
            input_dim = length(start_point);
            if (input_dim < 1) || (input_dim > 3)
                error('Input start_point must have 1 (in 1D), 2 (in 2D), or 3 (in 3D) elements.');
            end
            if input_dim ~= length(end_point)
                error('Inputs start_point and end_point must have the same number of elements.');
            end
            
            % check if this is the first element, and set the dimension
            if obj.number_elements == 0
                obj.dim = input_dim;
            end
            
            % check that the element is being added to an array with the
            % correct dimensions
            if obj.dim ~= input_dim
                error([num2str(input_dim) 'D line element cannot be added to an array with ' num2str(obj.dim) 'D elements.']);
            end           
            
            % increment the number of elements
            obj.number_elements = obj.number_elements + 1;
            
            % compute length
            line_length = norm(end_point - start_point);
            
            % add the element properties
            obj.elements{obj.number_elements}.group_id    = 0;
            obj.elements{obj.number_elements}.type        = 'line';
            obj.elements{obj.number_elements}.dim         = 1;
            obj.elements{obj.number_elements}.start_point = start_point(:).';
            obj.elements{obj.number_elements}.end_point   = end_point(:).';
            obj.elements{obj.number_elements}.active      = true;
            obj.elements{obj.number_elements}.measure     = line_length;
            
        end

        % remove element
        function removeElement(obj, element_num)
            
            % check the element number exists
            if element_num > obj.number_elements
                error(['Cannot remove element ' num2str(element_num) ' from array with ' num2str(obj.num_elements) ' elements']);
            end
            
            % remove the element
            obj.elements(element_num) = [];
            
            % decrement the number of elements
            obj.number_elements = obj.number_elements - 1;
            
        end
        
        % -----------------
        % WEIGHTS AND MASKS
        % -----------------
        
        % get grid weights for one element
        function grid_weights = getElementGridWeights(obj, kgrid, element_num)
            grid_weights = getOffGridPoints(obj, kgrid, element_num, false);
        end
        
        % get binary mask for one element
        function mask = getElementBinaryMask(obj, kgrid, element_num)
            mask = getOffGridPoints(obj, kgrid, element_num, true);
        end        
        
        % get grid weights for complete array
        function grid_weights = getArrayGridWeights(obj, kgrid)
            
            % check the array has elements
            obj.checkForElements(dbstack);
            
            % initialise the grid weights matrix
            grid_weights = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
            
            % loop through the elements and add element grid weights
            for ind = 1:obj.number_elements
                grid_weights = grid_weights + obj.getOffGridPoints(kgrid, ind, false);                
            end
            
        end          
        
        % get binary mask for complete array
        function mask = getArrayBinaryMask(obj, kgrid)
            
            % check the array has elements
            obj.checkForElements(dbstack);
            
            % initialise the mask
            mask = false(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
            
            % loop through the elements
            for ind = 1:obj.number_elements
            
                % get the offgrid source mask
                grid_weights = obj.getOffGridPoints(kgrid, ind, true);

                % add to binary mask
                mask = mask | grid_weights;
                
            end
            
        end        

        % compute distributed source signal
        function distributed_source_signal = getDistributedSourceSignal(obj, kgrid, source_signal)
            
            % update command line
            func_start_time = clock;
            disp('Computing distributed source signal...');
            
            % check the array has elements
            obj.checkForElements(dbstack);
            
            % get the binary mask to count how many points contribute to
            % the source
            comp_start_time = clock;
            fprintf('  calculating binary mask...                  ');
            mask = obj.getArrayBinaryMask(kgrid);
            mask_ind = find(mask);
            num_source_points = sum(mask(:));
            disp(['completed in ' scaleTime(etime(clock, comp_start_time))]);
          
            % number of time points in the signal
            Nt = size(source_signal, 2);
            
            % estimate size of the signal, and cast if needed
            if obj.single_precision
                data_type = 'single';
                sz_bytes = num_source_points * Nt * 4;
            else
                data_type = 'double';
                sz_bytes = num_source_points * Nt * 8;
            end
            sz_ind = 1;
            while sz_bytes > 1024
                sz_bytes = sz_bytes / 1024;
                sz_ind = sz_ind + 1;
            end
            prefixes = {'', 'K', 'M', 'G', 'T'};
            sz_bytes = round(sz_bytes, 2, 'significant');
            disp(['  approximate size of source matrix:          ' num2str(sz_bytes) ' ' prefixes{sz_ind} 'B (' data_type ' precision)']);
            
            % cast input source signal
            source_signal = cast(source_signal, data_type);
            
            % initialise the source signal
            distributed_source_signal = zeros(num_source_points, Nt, data_type);
            
            % loop through the elements
            for ind = 1:obj.number_elements
            
                % get the offgrid source weights
                comp_start_time = clock;
                fprintf(['  calculating element ' num2str(ind) ' grid weights...       ']);
                source_weights = obj.getElementGridWeights(kgrid, ind);
                disp(['completed in ' scaleTime(etime(clock, comp_start_time))]);
                
                % get indices of the non-zero points 
                element_mask_ind = find(source_weights ~= 0);
                
                % convert these to indices in the distributed source
                local_ind = ismember(mask_ind, element_mask_ind);
                
                % add to distributed source
                comp_start_time = clock;
                fprintf(['  calculating element ' num2str(ind) ' distributed source... ']);
                distributed_source_signal(local_ind, :) = ...
                    distributed_source_signal(local_ind, :) ...
                    + bsxfun(@times, source_weights(element_mask_ind), source_signal(ind, :));
                disp(['completed in ' scaleTime(etime(clock, comp_start_time))]);
                
            end
            
            disp(['  total computation time ' scaleTime(etime(clock, func_start_time))]);
            
        end
        
        % combine distributed sensor data
        function combined_sensor_data = combineSensorData(obj, kgrid, sensor_data)
            
            % check the array has elements
            obj.checkForElements(dbstack);
            
            % get the binary mask and the indices of the active points
            mask = obj.getArrayBinaryMask(kgrid);
            mask_ind = find(mask);
            
            % number of time points in the signal
            Nt = size(sensor_data, 2);
            
            % initialise the combined sensor data
            combined_sensor_data = zeros(obj.number_elements, Nt);
            
            % loop through the array elements
            for element_num = 1:obj.number_elements
            
                % get the offgrid source weights
                source_weights = obj.getElementGridWeights(kgrid, element_num);
                
                % get indices of the non-zero points 
                element_mask_ind = find(source_weights ~= 0);
                
                % convert these to indices in the sensor data
                local_ind = ismember(mask_ind, element_mask_ind);
                
                % get the sensor data that belongs to this element, weight
                % by the source weights, and sum 
                combined_sensor_data(element_num, :) = sum(bsxfun(@times, sensor_data(local_ind, :), source_weights(element_mask_ind)), 1);
                
                % compute measure (length/area/volume) in grid squares (assuming dx = dy = dz)
                m_grid = obj.elements{element_num}.measure ./ (kgrid.dx).^(obj.elements{element_num}.dim);
                
                % normalise by the element measure
                combined_sensor_data(element_num, :) = combined_sensor_data(element_num, :) ./ m_grid;
                
                % this will need some additional weighting factor to give
                % the sensor data in pressure
                
                % consider what to do if the elements are not active?
                                
            end
            
        end
        
        % set the position and orientation of the array
        function setArrayPosition(obj, translation, rotation)
            
            % get affine matrix and store
            obj.array_transformation = getAffineMatrix(translation, rotation);
            
        end
        
        % set the affine transform to specify the position and orientation
        % of the array
        function setAffineTransform(obj, affine_transform)
            
            % check array has elements
            if obj.number_elements == 0
                error('Array must have at least one element before a transformation can be defined.');
            end
            
            % check the input is the correct size
            [sx, sy] = size(affine_transform);
            if (obj.dim == 2 && (sx ~= 3 || sy ~= 3)) || (obj.dim == 3 && (sx ~= 4 || sy ~= 4))
                error('Affine transform must be a 3x3 matrix for arrays in 2D, and 4x4 for arrays in 3D.'); 
            end
            
            % assign the transform
            obj.array_transformation = affine_transform;
            
        end
            
        % get the positions of the elements after transformation
        function element_pos = getElementPositions(obj)
            
            % initialise output
            element_pos = zeros(obj.dim, obj.number_elements);
            
            % loop through the elements and assign transformed position
            for element_num = 1:obj.number_elements
                element_pos(:, element_num) = obj.affine(obj.elements{element_num}.position);
            end
            
        end
        
        % set the optional input parameters
        function setOptionalInputs(obj, varargin)
    
            % check inputs are given as pairs
            if rem(length(varargin), 2)
                error('Optional input parameters must be given as param, value pairs.');
            end
            
            % loop through the optional inputs
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'Axisymmetric'
                            
                            % assign input
                            obj.axisymmetric = varargin{input_index + 1};
                            
                        case 'BLITolerance'

                            % assign input
                            obj.bli_tolerance = varargin{input_index + 1};

                            % check parameter range
                            if ~isnumeric(obj.bli_tolerance) || (obj.bli_tolerance < 0) || (obj.bli_tolerance > 1)
                                error('Optional input ''BLITolerance'' must be between 0 and 1.');
                            end

                        case 'BLIType'

                            % assign input
                            obj.bli_type = varargin{input_index + 1};

                            % check parameter
                            if ~any(strcmp(obj.bli_type, {'sinc', 'exact'}))
                                error('Optional input ''BLIType'' must be either ''sinc'' or ''exact''.')
                            end

                        case 'SinglePrecision'
                            
                            % assign input
                            obj.single_precision = varargin{input_index + 1};

                            % check parameter
                            if ~islogical(obj.single_precision)
                                error('Optional input ''SinglePrecision'' must be Boolean.');
                            end  
                            
                        case 'UpsamplingRate'

                            % assign input
                            obj.upsampling_rate = varargin{input_index + 1};

                        otherwise
                            error('Unknown optional input.');
                    end
                end
            end
            
        end
        
        % plot the array
        function plotArray(obj, new_figure)

            % create new figure if needed
            if (nargin == 1) || isempty(new_figure) || new_figure
                figure;
            end
            hold on;
            
            % loop through the elements
            for element_num = 1:obj.number_elements
                switch obj.elements{element_num}.type
                    case 'bowl'

                    case 'disc'
                        
                        % Disc plotting code adapted from MATHWORKS file exchange submission
                        % https://uk.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d
                        %
                        % Copyright (c) 2010, Christian Reinbacher
                        % All rights reserved.
                        %
                        % Redistribution and use in source and binary forms, with or without
                        % modification, are permitted provided that the following conditions are met:
                        %
                        % * Redistributions of source code must retain the above copyright notice, this
                        %   list of conditions and the following disclaimer.
                        %
                        % * Redistributions in binary form must reproduce the above copyright notice,
                        %   this list of conditions and the following disclaimer in the documentation
                        %   and/or other materials provided with the distribution
                        %
                        % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
                        % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
                        % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
                        % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
                        % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
                        % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
                        % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
                        % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
                        % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
                        % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                        
                        % get points on the perimeter of the disc
                        angles = 0:0.05:2*pi;
                        v = null(obj.affine(obj.elements{element_num}.focus_position) - obj.affine(obj.elements{element_num}.position));
                        points = repmat(obj.affine(obj.elements{element_num}.position).', 1, length(angles)) + obj.elements{element_num}.diameter/2 * (v(:, 1) * cos(angles) + v(:, 2) * sin(angles));
                        
                        % create filled disc
                        fill3(points(1, :), points(2, :), points(3, :), obj.element_plot_colour);
                        
                    case 'arc'
                        
                        % get arc
                        arc = makeCartArc(...,
                            obj.affine(obj.elements{element_num}.position), ...
                            obj.elements{element_num}.radius_of_curvature, ...
                            obj.elements{element_num}.diameter, ...
                            obj.affine(obj.elements{element_num}.focus_position), ...
                            obj.num_arc_plot_points);
                        
                        % plot
                        plot(arc(2, :), arc(1, :), '-', 'Color', obj.element_plot_colour);
                        
                    otherwise
                        error([obj.elements{element_num}.type ' is not a valid array element type.']);
                end
            end
            
            % annotate plot
            switch obj.dim
                case 2
                    axis image;
                    axis ij;
                    box on;
                    xlabel('y [m]');
                    ylabel('x [m]');
                case 3
                    axis image;
                    xlabel('x [m]');
                    ylabel('y [m]');
                    zlabel('z [m]');
            end
        end
            
    end
    
    % internal class methods only accessible by other functions 
    methods (Hidden = true, Access = 'protected') 
        
        % check if the array has elements
        function checkForElements(obj, call_stack)
            if obj.number_elements == 0
                error(['Cannot call ' call_stack.name ' on an array with zero elements.']);
            end
        end
        
        % apply affine transformation
        function vec = affine(obj, vec)
            
            if isempty(obj.array_transformation)
                return;
            end
               
            % check vector is the correct length
            if obj.dim ~= numel(vec)
                error('Input vector length must match the number of dimensions');
            end
            
            % apply transformation
            vec = obj.array_transformation * [vec(:); 1];
            vec = vec(1:end-1)';
            
        end
        
        % get grid weights or mask for one element
        function grid_weights = getOffGridPoints(obj, kgrid, element_num, mask_only)

            % check the array has elements
            obj.checkForElements(dbstack);
            
            % check inputs
            validateattributes(kgrid, {'kWaveGrid'}, {}, 'getElementGridWeights', 'kgrid', 2);
            validateattributes(element_num, {'numeric'}, {'integer', '>=', 1, '<=', obj.number_elements}, 'getElementGridWeights', 'element_num', 2);
            
            % compute measure (length/area/volume) in grid squares (assuming dx = dy = dz)
            m_grid = obj.elements{element_num}.measure ./ (kgrid.dx).^(obj.elements{element_num}.dim);
            
            % get number of integration points
            if strcmp(obj.elements{element_num}.type, 'custom')
                
                % assign number of integration points directly
                m_integration = size(obj.elements{element_num}.integration_points, 2);
                
            else
                
                % compute the number of integration points using the upsampling rate          
                m_integration = ceil(m_grid .* obj.upsampling_rate);
                
            end
            
            % compute integration points covering element
            switch obj.elements{element_num}.type
                case 'annulus'
                    
                    % compute points using makeCartSphericalSegment
                    integration_points = makeCartSphericalSegment(...
                        obj.affine(obj.elements{element_num}.position), ...
                        obj.elements{element_num}.radius_of_curvature, ...
                        obj.elements{element_num}.inner_diameter, ...
                        obj.elements{element_num}.outer_diameter, ...
                        obj.affine(obj.elements{element_num}.focus_position), ...
                        m_integration);
                    
                case 'arc'
                    
                    % compute points using makeCartArc
                    integration_points = makeCartArc(...
                        obj.affine(obj.elements{element_num}.position), ...
                        obj.elements{element_num}.radius_of_curvature, ...
                        obj.elements{element_num}.diameter, ...
                        obj.affine(obj.elements{element_num}.focus_position), ...
                        m_integration);
                    
                case 'bowl'
                    
                    % compute points using makeCartBowl
                    integration_points = makeCartBowl(...
                        obj.affine(obj.elements{element_num}.position), ...
                        obj.elements{element_num}.radius_of_curvature, ...
                        obj.elements{element_num}.diameter, ...
                        obj.affine(obj.elements{element_num}.focus_position), ...
                        m_integration);
                    
                case 'custom'
                    
                    % directly assign integration points
                    integration_points = obj.elements{element_num}.integration_points;
                    
                case 'disc'
                    
                    % compute points using makeCartDisc
                    integration_points = makeCartDisc(...
                        obj.affine(obj.elements{element_num}.position), ...
                        obj.elements{element_num}.diameter/2, ...
                        obj.affine(obj.elements{element_num}.focus_position), ...
                        m_integration, ...
                        false, ...
                        obj.use_spiral_disc_points);
                    
                case 'rect'
                    
                    % compute points using makeCartRect
                    integration_points = makeCartRect(...
                        obj.affine(obj.elements{element_num}.position), ...
                        obj.elements{element_num}.length, ...
                        obj.elements{element_num}.width, ...
                        obj.elements{element_num}.orientation, ...
                        m_integration);
                    
                case 'line'
                    
                    % get distance between points in each dimension
                    d = (obj.elements{element_num}.end_point - obj.elements{element_num}.start_point) ./ m_integration;

                    % compute a set of uniformly spaced Cartesian points
                    % covering the line using linspace, where the end
                    % points are offset by half the point spacing
                    switch obj.dim
                        case 1
                            integration_points = linspace(obj.elements{element_num}.start_point + d(1)/2, obj.elements{element_num}.end_point - d(1)/2, m_integration);
                        case 2
                            px = linspace(obj.elements{element_num}.start_point(1) + d(1)/2, obj.elements{element_num}.end_point(1) - d(1)/2, m_integration);
                            py = linspace(obj.elements{element_num}.start_point(2) + d(2)/2, obj.elements{element_num}.end_point(2) - d(2)/2, m_integration);
                            integration_points = [px; py];
                        case 3
                            px = linspace(obj.elements{element_num}.start_point(1) + d(1)/2, obj.elements{element_num}.end_point(1) - d(1)/2, m_integration);
                            py = linspace(obj.elements{element_num}.start_point(2) + d(2)/2, obj.elements{element_num}.end_point(2) - d(2)/2, m_integration);
                            pz = linspace(obj.elements{element_num}.start_point(3) + d(3)/2, obj.elements{element_num}.end_point(3) - d(3)/2, m_integration);
                            integration_points = [px; py; pz];
                    end
                    
                otherwise
                    
                    % unknown element type
                    error([obj.elements{element_num}.type ' is not a valid array element type.']);
                    
            end
            
            % recompute actual number of points
            m_integration = size(integration_points, 2);
            
            % compute scaling factor
            scale = m_grid ./ m_integration;
            
            if obj.axisymmetric
            
                % create new expanded grid
                kgrid_expanded = kWaveGrid(kgrid.Nx, kgrid.dx, 2 * kgrid.Ny, kgrid.dy);
                
                % remove integration points which are outside grid
                integration_points = trimCartPoints(kgrid_expanded, integration_points);
                
                % calculate grid weights from BLIs centered on the integration points
                grid_weights = offGridPoints(kgrid_expanded, integration_points, scale, ...
                    'BLITolerance', obj.bli_tolerance, ...
                    'BLIType', obj.bli_type, ...
                    'MaskOnly', mask_only, ...
                    'SinglePrecision', obj.single_precision);
                
                % keep points in the positive y domain
                grid_weights  = grid_weights(:, kgrid.Ny + 1:end);
                
            else
                
                % remove integration points which are outside grid
                integration_points = trimCartPoints(kgrid, integration_points);
                
                % calculate grid weights from BLIs centered on the integration points
                grid_weights = offGridPoints(kgrid, integration_points, scale, ...
                    'BLITolerance', obj.bli_tolerance, ...
                    'BLIType', obj.bli_type, ...
                    'MaskOnly', mask_only, ...
                    'SinglePrecision', obj.single_precision);
                
            end
            
        end  
        
    end
    
end