%KWAVEGRID Class definition for k-Wave grid.
%
% DESCRIPTION:
%     kWaveGrid is the grid class used across the k-Wave Toolbox. An object
%     of the kWaveGrid class contains the grid coordinates and wavenumber
%     matrices used within the simulation and reconstruction functions in
%     k-Wave. The grid matrices are indexed as: (x, 1) in 1D; (x, y) in
%     2D; and (x, y, z) in 3D. The grid is assumed to be a regularly spaced
%     Cartesian grid, with grid spacing given by dx, dy, dz (typically the
%     grid spacing in each direction is constant). 
%
% USAGE:
%     kgrid = kWaveGrid(Nx, dx)
%     kgrid = kWaveGrid(Nx, dx, Ny, dy)
%     kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz)
%
% INPUTS:
%     Nx, Ny, Nz      - number of grid points in each Cartesian direction
%     dx, dy, dz      - grid point spacing in each Cartesian direction [m]
%
% OUTPUTS:
%     kgrid           - kWaveGrid object used by the simulation and 
%                       reconstructions functions within k-Wave
%
% DYNAMIC PROPERTIES:
%     Properties which can be queried or modified after the object is
%     created. 
%
%     .Nt             - number of time steps (default = 'auto')
%     .dt             - time step [s] (default = 'auto')
%     .t_array        - evenly spaced array of time values [s] (default =
%                       'auto')
%
%     Note, t_array is a derived property, thus changing Nt or dt will
%     also change t_array, and modifying t_array will update Nt and dt.
%
% STATIC PROPERTIES:
%     Properties which can be queried, but not modified, after the object
%     is created.  
%
%     .k              - Nx x Ny x Nz grid of the scalar wavenumber, where
%                       k = sqrt(kx.^2 + ky.^2 + kz.^2) [rad/m]
%     .k_max          - maximum wavenumber (spatial frequency) supported
%                       by the grid [rad/m] 
%     .dim            - number of spatial dimensions
%     .total_grid_points 
%                     - total number of grid points (equal to Nx * Ny * Nz)
%     .highest_prime_factors 
%                     - 3 x 1 vector of the highest prime factor in each
%                       direction  
%
%     And for each spatial dimension x, y, z:
%
%     .Nx             - number of grid points in the x-direction
%     .dx             - grid point spacing in the x-direction [m]
%     .x_vec          - Nx x 1 vector of the grid coordinates in the
%                       x-direction [m] 
%     .x              - Nx x Ny x Nz grid containing repeated copies of the
%                       grid coordinates in the x-direction [m] 
%     .x_size         - length of grid dimension in the x-direction [m]
%     .kx_vec         - Nx x 1 vector of the wavenumber components in the
%                       x-direction [rad/m]
%     .kx             - Nx x Ny x Nz grid containing repeated copies of the
%                       wavenumber components in the x-direction [rad/m]
%     .kx_max         - maximum wavenumber (spatial frequency) supported by
%                       the grid in the x-direction [rad/m] 
%
% METHODS - SETTIME:
%     setTime sets dt and Nt to the input values. The syntax is:
%
%         .setTime(Nt, dt)  
%
% METHODS - MAKETIME:
%     makeTime automatically specifies Nt, dt, and t_array based on the
%     Courant-Friedrichs-Lewy (CFL) number and the grid size. The sytax is:
%
%         .makeTime(sound_speed)
%         .makeTime(sound_speed, cfl)
%         .makeTime(sound_speed, cfl, t_end)
%         .makeTime(sound_speed, [], t_end)
%
%     where
%
%         sound_speed - sound speed in the medium [m/s]
%         cfl         - CFL number (default = 0.3)
%         t_end       - simulation time [s]
%
%     The time-step dt is chosen based on the CFL number (the default is
%     0.3), where dt = cfl * dx / sound_speed. If t_end is not specified,
%     the number of time-steps is chosen based on the time it takes to
%     travel from one corner of the grid to the geometrically opposite
%     corner. Note, if sound_speed is given as a matrix, the calculation
%     for dt is based on the maximum value, and the calculation for t_end
%     based on the minimum value. 
%
% METHODS - DISCRETE TRIGONOMETRIC TRANSFORMS:
%     Objects of the kWaveGrid class also have the following methods to
%     return the wavenumbers for use with discrete sine and cosine
%     transforms. The syntax is:
%
%         k      = kgrid.k_dtt(dtt_type)
%         kx_vec = kgrid.kx_vec_dtt(dtt_type)
%         ky_vec = kgrid.ky_vec_dtt(dtt_type)
%         kz_vec = kgrid.kz_vec_dtt(dtt_type)
%
%     where
%
%         k           - Nx x Ny x Nz grid of the scalar wavenumber, where 
%                       k = sqrt(kx.^2 + ky.^2 + kz.^2) [rad/m] 
%         kx_vec      - Nx x 1 vector of the wavenumber components in the
%                       x-direction [rad/m] 
%         ky_vec      - Ny x 1 vector of the wavenumber components in the
%                       y-direction [rad/m] 
%         kz_vec      - Nz x 1 vector of the wavenumber components in the
%                       z-direction [rad/m] 
%
%     Here dtt_type is the type of discrete trigonometric transform (this
%     determines the assumed symmetry of the input function), where:  
%
%         1 -> DCT-I    WSWS (W: Whole Sample, S: Symmetric)
%         2 -> DCT-II   HSHS 
%         3 -> DCT-III  WSWA
%         4 -> DCT-IV   HSHA
%         5 -> DST-I    WAWA
%         6 -> DST-II   HAHA (H: Half Sample, A: Antisymmetric)
%         7 -> DST-III  WAWS
%         8 -> DST-IV   HAHS
%
%     These methods can also return the implied periodic length of the
%     grid (denoted M) given the grid size and the dtt_type, e.g.,
%
%         [kx_vec, Mx] = kgrid.kx_vec_dtt(1);
%         [k, M]       = kgrid.k_dtt(1);
%
%     The method .k_dtt returns the product of the implied length in each
%     dimension, e.g., in 2D, M = Mx*My. This is used for scaling the
%     inverse trigonometric transforms.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 22nd July 2010
%     last update - 8th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby
%
% See also kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kspaceFirstOrderAS, kspaceSecondOrder, pstdElastic2D, pstdElastic2D

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
classdef kWaveGrid < handle
    
    % define literals
    properties (Constant, Hidden = true)
        
        % default CFL number
        CFL_DEFAULT = 0.3;
        
        % machine precision
        MACHINE_PRECISION = 100 * eps;
        
    end
    
    % define the private properties (these parameters are stored and cannot
    % be altered by the user after the object is created, but they can be
    % accessed). The numbers assigned here are the default values used if
    % the parameters are not set by the constructor.
    properties (GetAccess = public, SetAccess = private)
        
        % grid size in x-direction [grid points]
        Nx = 0;
        
        % grid size in y-direction [grid points]
        Ny = 0;
        
        % grid size in z-direction [grid points]
        Nz = 0;
        
        % grid point spacing in x-direction [m]
      	dx = 0;
        
        % grid point spacing in y-direction [m]
        dy = 0;
        
        % grid point spacing in z-direction [m]
        dz = 0;
        
        % Nx x 1 vector of wavenumber components in the x-direction [rad/m]
        kx_vec = 0;
        
        % Ny x 1 vector of wavenumber components in the y-direction [rad/m]
        ky_vec = 0;
        
        % Nz x 1 vector of wavenumber components in the z-direction [rad/m]
        kz_vec = 0;           
        
        % scalar wavenumber [rad/m]
        k = 0;
        
        % maximum supported spatial frequency in the x-direction [rad/m]
        kx_max = 0;
        
        % maximum supported spatial frequency in the y-direction [rad/m]
        ky_max = 0;
        
        % maximum supported spatial frequency in the z-direction [rad/m]
        kz_max = 0;
        
        % maximum supported spatial frequency in all directions [rad/m]
        k_max = 0;
        
        % number of dimensions
        dim = 0;      
        
        % non-uniform grid
        nonuniform = false;
        
    end
    
    % define the public properties (these parameters are stored and can be
    % both accessed and altered by the user after the object is created)
    properties (GetAccess = public, SetAccess = public)
        
        % size of time step [s]
        dt = 'auto';
        
        % number of time steps [s]
        Nt = 'auto';
                
    end    
    
    % define the dependent properties without set methods (these parameters
    % are not stored but re-computed each time they are needed)
    properties(Dependent = true, GetAccess = public, SetAccess = private)
        
        % Nx x 1 vector of the grid coordinates in the x-direction [m] 
        x_vec;
        
        % Ny x 1 vector of the grid coordinates in the y-direction [m] 
        y_vec;
        
        % Nz x 1 vector of the grid coordinates in the z-direction [m] 
        z_vec;    
        
        % Nx x Ny x Nz grid containing repeated copies of the grid coordinates in the x-direction [m]  
        x;
        
        % Nx x Ny x Nz grid containing repeated copies of the grid coordinates in the y-direction [m]  
        y;
        
        % Nx x Ny x Nz grid containing repeated copies of the grid coordinates in the z-direction [m]  
        z;
        
        % Nx x Ny x Nz grid containing repeated copies of the wavenumber components in the x-direction [rad/m] 
        kx;
        
        % Nx x Ny x Nz grid containing repeated copies of the wavenumber components in the y-direction [rad/m] 
        ky;
        
        % Nx x Ny x Nz grid containing repeated copies of the wavenumber components in the z-direction [rad/m] 
        kz;
        
        % size of grid in the x-direction [m]
        x_size;
        
        % size of grid in the y-direction [m]
        y_size;
        
        % size of grid in the z-direction [m]
        z_size;
        
        % total number of grid points
        total_grid_points;
        
    end
    
    % define the dependent properties with a set method (these parameters
    % are not stored but re-computed each time they are needed, and the set
    % method modifies other properties) 
    properties(Dependent = true, GetAccess = public, SetAccess = public)
    
        % time array [s]
        t_array;
        
    end    
    
    % define the staggered grid properties (these can be accessed but not
    % altered by the user after the object is created, the parameters are
    % also hidden from view)
    properties (Hidden = true, GetAccess = public, SetAccess = private)
        
        % position vectors for the staggered grid points in [0, 1]
        xn_vec = 0;
        yn_vec = 0;
        zn_vec = 0;
        xn_vec_sgx = 0;
        yn_vec_sgy = 0;
        zn_vec_sgz = 0;    
        
        % transformation gradients between uniform and staggered grids
        dxudxn = 0;
        dyudyn = 0;
        dzudzn = 0;
        dxudxn_sgx = 0;
        dyudyn_sgy = 0;
        dzudzn_sgz = 0;  
        
    end
    
    % define the dependent staggered grid properties (these parameters are
    % not stored but re-computed each time they are needed, the parameters
    % are also hidden from view)
    properties(Dependent = true, Hidden = true, GetAccess = public, SetAccess = private)
        
        % 3D plaid non-uniform spatial grids
        xn;
        yn;
        zn;
        
    end
    
    % constructor function
    methods
        function kgrid = kWaveGrid(varargin)
                        
            % assign the input values to the grid object
            switch nargin
                case 6
                        
                    % 3D uniform grid
                    kgrid.Nx = varargin{1};
                    kgrid.dx = varargin{2};
                    kgrid.Ny = varargin{3};
                    kgrid.dy = varargin{4};
                    kgrid.Nz = varargin{5};
                    kgrid.dz = varargin{6};

                    % set the number of dimensions
                    kgrid.dim = 3;
            
                case 4
                    
                    % 2D uniform grid
                    kgrid.Nx = varargin{1};
                    kgrid.dx = varargin{2};
                    kgrid.Ny = varargin{3};
                    kgrid.dy = varargin{4};

                    % set the number of dimensions
                    kgrid.dim = 2;
                    
                case 2
                    
                    % 1D uniform grid
                    kgrid.Nx = varargin{1};
                    kgrid.dx = varargin{2};
                    
                    % set the number of dimensions
                    kgrid.dim = 1;
                    
                otherwise
                    error('Incorrect number of input arguments');
            end
                       
            switch kgrid.dim
                case 1
                    
                    % assign the grid parameters for the x spatial direction
                    kgrid.kx_vec = kgrid.makeDim(kgrid.Nx, kgrid.dx);
                   
                    % define the scalar wavenumber based on the wavenumber components
                    kgrid.k = abs(kgrid.kx_vec);
                    
                    % define maximum supported frequency
                    kgrid.kx_max = max(abs(kgrid.kx_vec(:)));
                    kgrid.k_max  = kgrid.kx_max;
                    
                case 2
                    
                    % assign the grid parameters for the x and y spatial directions
                    kgrid.kx_vec = kgrid.makeDim(kgrid.Nx, kgrid.dx);
                    kgrid.ky_vec = kgrid.makeDim(kgrid.Ny, kgrid.dy);

                    % define the wavenumber based on the wavenumber components
                    % (using bsxfun saves memory by avoiding ndgrid)
                    kgrid.k = zeros(kgrid.Nx, kgrid.Ny);
                    kgrid.k = bsxfun(@plus, (reshape(kgrid.kx_vec, [], 1, 1).^2), kgrid.k);
                    kgrid.k = bsxfun(@plus, (reshape(kgrid.ky_vec, 1, [], 1).^2), kgrid.k);
                    kgrid.k = sqrt(kgrid.k);
                    
                    % define maximum supported frequency
                    kgrid.kx_max = max(abs(kgrid.kx_vec(:)));
                    kgrid.ky_max = max(abs(kgrid.ky_vec(:)));
                    kgrid.k_max  = min([kgrid.kx_max, kgrid.ky_max]); 
                    
                case 3
                    
                    % assign the grid parameters for the x, y, and z spatial directions
                    kgrid.kx_vec = kgrid.makeDim(kgrid.Nx, kgrid.dx);
                    kgrid.ky_vec = kgrid.makeDim(kgrid.Ny, kgrid.dy);
                    kgrid.kz_vec = kgrid.makeDim(kgrid.Nz, kgrid.dz);

                    % define the wavenumber based on the wavenumber components
                    % (using bsxfun saves memory by avoiding ndgrid)
                    kgrid.k = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
                    kgrid.k = bsxfun(@plus, (reshape(kgrid.kx_vec, [], 1, 1).^2), kgrid.k);
                    kgrid.k = bsxfun(@plus, (reshape(kgrid.ky_vec, 1, [], 1).^2), kgrid.k);
                    kgrid.k = bsxfun(@plus, (reshape(kgrid.kz_vec, 1, 1, []).^2), kgrid.k);
                    kgrid.k = sqrt(kgrid.k);                   
                   
                    % define maximum supported frequency
                    kgrid.kx_max = max(abs(kgrid.kx_vec(:)));
                    kgrid.ky_max = max(abs(kgrid.ky_vec(:)));
                    kgrid.kz_max = max(abs(kgrid.kz_vec(:)));
                    kgrid.k_max  = min([kgrid.kx_max, kgrid.ky_max, kgrid.kz_max]);  
                    
            end
        end
    end
    
    % functions for dependent variables that only run when queried
    methods
        
        % calculate x based on kx
        function x = get.x(obj)
            x = obj.x_size .* obj.kx .* obj.dx ./ (2 .* pi);
        end
        
        % calculate y based on ky
        function y = get.y(obj)
            y = obj.y_size .* obj.ky .* obj.dy ./ (2 .* pi);
        end
        
        % calculate z based on kz
        function z = get.z(obj)
            z = obj.z_size .* obj.kz .* obj.dz ./ (2 .* pi);
        end        
        
        % duplicate kx vector using repmat
        function kx = get.kx(obj)
            switch obj.dim
                case 1
                    kx = obj.kx_vec;
                case 2
                    kx = repmat(obj.kx_vec, [1, obj.Ny]);
                case 3
                    kx = repmat(obj.kx_vec, [1, obj.Ny, obj.Nz]);
            end
        end
        
        % rotate ky vector to y direction then duplicate using repmat
        function ky = get.ky(obj)
            switch obj.dim
                case 1
                    ky = 0;
                case 2
                    ky = repmat(obj.ky_vec.', [obj.Nx, 1]);
                case 3
                    ky = repmat(obj.ky_vec.', [obj.Nx, 1, obj.Nz]);
            end
        end     
        
        % permute kz vector to z direction then duplicate using repmat
        function kz = get.kz(obj)
            switch obj.dim
                case 1
                    kz = 0;
                case 2
                    kz = 0;
                case 3
                    kz = repmat(permute(obj.kz_vec, [2 3 1]), [obj.Nx, obj.Ny, 1]);
            end
        end
  
        % calculate x_vec based on kx_vec
        function x_vec = get.x_vec(obj)
            x_vec = obj.x_size .* obj.kx_vec .* obj.dx ./ (2 .* pi);
        end
        
        % calculate y_vec based on ky_vec
        function y_vec = get.y_vec(obj)
            y_vec = obj.y_size .* obj.ky_vec .* obj.dy ./ (2 .* pi);
        end
        
        % calculate z_vec based on kz_vec
        function z_vec = get.z_vec(obj)
            z_vec = obj.z_size .* obj.kz_vec .* obj.dz ./ (2 .* pi);
        end        
                        
        % calculate x_size based on Nx and dx
        function x_size = get.x_size(obj)
            x_size = obj.Nx .* obj.dx;
        end
        
        % calculate y_size based on Ny and dy
        function y_size = get.y_size(obj)
            y_size = obj.Ny .* obj.dy;
        end
        
        % calculate z_size based on Nz and dz
        function z_size = get.z_size(obj)
            z_size = obj.Nz .* obj.dz;
        end     
        
        % calulate t_array based on Nt and dt
        function t_array = get.t_array(obj)
            if strcmp(obj.Nt, 'auto') || strcmp(obj.dt, 'auto')
                t_array = 'auto';
            else
                t_array = (0:obj.Nt - 1) .* obj.dt;
            end
        end
        
        % set Nt and dt based on t_array
        function set.t_array(obj, t_array)
            
            % check for 'auto' input
            if strcmp(t_array, 'auto')
                
                % set values to auto
                obj.Nt = 'auto';
                obj.dt = 'auto';
                
            else

                % extract property values
                Nt_temp = length(t_array);
                dt_temp = t_array(2) - t_array(1);

                % check the time array begins at zero
                if t_array(1) ~= 0
                    error('t_array must begin at zero.');
                end
                
                % check the time array is evenly spaced
                if max(t_array(2:end) - t_array(1:end-1) - dt_temp) > obj.MACHINE_PRECISION
                    error('t_array must be evenly spaced.');
                end
                
                % check the time steps are increasing
                if dt_temp <= 0
                    error('t_array must be monotonically increasing.');
                end
                
                % assign values
                obj.Nt = Nt_temp;
                obj.dt = dt_temp;                

            end
        end
        
        % calculate the total number of points in the grid
        function N = get.total_grid_points(obj)
            switch obj.dim
                case 1
                    N = obj.Nx;
                case 2
                    N = obj.Nx .* obj.Ny;
                case 3
                    N = obj.Nx .* obj.Ny .* obj.Nz;
            end
        end
        
    end
       
    % general class methods
    methods
       
        % set Nt and dt based on user input
        function setTime(obj, Nt, dt)
            
            % check the value for Nt
            if (round(Nt) ~= Nt) || (Nt <= 0)
                error('Nt must be a positive integer.');
            end
            
            % check the value for dt
            if dt <= 0
                error('dt must be positive.');
            end
            
            % assign values
            obj.Nt = Nt;
            obj.dt = dt;
            
        end
        
        % compute Nt and dt based on the cfl number and grid size, where
        % the number of time-steps is chosen based on the time it takes to
        % travel from one corner of the grid to the geometrically opposite
        % corner. Note, if c is given as a matrix, the calculation for dt
        % is based on the maximum value, and the calculation for t_end
        % based on the minimum value. 
        function [t_array, dt] = makeTime(obj, c, cfl, t_end)
            
            % check user defined value for the Courant-Friedrichs-Lewy
            % stability number , otherwise assign default value
            if nargin < 3 || isempty(cfl)
                cfl = obj.CFL_DEFAULT;
            end

            % if c is a matrix, find the minimum and maximum values
            c_max = max(c(:));
            c_min = min(c(:));

            % check for user define t_end, otherwise set the simulation
            % length based on the size of the grid diagonal and the maximum
            % sound speed in the medium 
            if nargin < 4
                switch obj.dim
                    case 1
                        t_end = obj.x_size ./ c_min;
                    case 2
                        t_end = sqrt(obj.x_size.^2 + obj.y_size.^2) ./ c_min;
                    case 3
                        t_end = sqrt(obj.x_size.^2 + obj.y_size.^2 + obj.z_size.^2) ./ c_min;
                end
            end

            % extract the smallest grid spacing
            switch obj.dim
                case 1
                    min_grid_dim = obj.dx;
                case 2
                    min_grid_dim = min([obj.dx, obj.dy]);
                case 3
                    min_grid_dim = min([obj.dx, obj.dy, obj.dz]);
            end

            % assign time step based on CFL stability criterion
            obj.dt = cfl .* min_grid_dim ./ c_max;
            
            % assign number of time steps based on t_end
            obj.Nt = floor(t_end / obj.dt) + 1;

            % catch case were dt is a recurring number
            if (floor(t_end / obj.dt) ~= ceil(t_end / obj.dt)) && (rem(t_end, obj.dt) == 0)
               obj.Nt = obj.Nt + 1; 
            end
            
            % assign outputs when using function as (k-Wave <= v1.1):
            %
            %     kgrid.t_array = makeTime(kgrid, c, ...);
            %
            % note, this produces a circular assignment, where the makeTime
            % method creates t_array, which is then copied to the t_array
            % of the same object, however, this allows backwards
            % compatability
            if nargout == 2
                t_array = obj.t_array;
                dt = obj.dt;
            elseif nargout == 1
                t_array = obj.t_array;
            end
            
        end
        
        % compute the individual wavenumber vectors, where dtt_type is the
        % type of discrete trigonometric transform, which corresponds to
        % the assumed input symmetry of the input function, where:  
        %    1: DCT-I    WSWS
        %    2: DCT-II   HSHS
        %    3: DCT-III  WSWA
        %    4: DCT-IV   HSHA
        %    5: DST-I    WAWA
        %    6: DST-II   HAHA
        %    7: DST-III  WAWS
        %    8: DST-IV   HAHS
        function [k, M] = k_dtt(obj, dtt_type)
           
            % check dtt_type is a scalar or a vector the same size as
            % obj.dim
            if ~( numel(dtt_type) ==  1 || numel(dtt_type) ==  obj.dim )
                error(['dtt_type must be a scalar, or ' num2str(obj.dim) 'D vector']);
            end
            
            switch obj.dim
                case 1
                    
                    % assign the grid parameters for the x direction
                    [k, M] = obj.kx_vec_dtt(dtt_type);
                    
                case 2
                    
                    % assign the grid parameters for the x and y spatial
                    % directions  
                    [kx_vec_dtt, Mx] = obj.kx_vec_dtt(dtt_type(1));
                    [ky_vec_dtt, My] = obj.ky_vec_dtt(dtt_type(end));

                    % define the wavenumber based on the wavenumber
                    % components
                    k = zeros(obj.Nx, obj.Ny);
                    k = bsxfun(@plus, (reshape(kx_vec_dtt, [], 1, 1).^2), k);
                    k = bsxfun(@plus, (reshape(ky_vec_dtt, 1, [], 1).^2), k);
                    k = sqrt(k);
                    
                    % define product of implied period
                    M = Mx*My;
                    
                case 3
                    
                    % assign the grid parameters for the x, y, and z
                    % spatial directions 
                    [kx_vec_dtt, Mx] = obj.kx_vec_dtt(dtt_type(1));
                    [ky_vec_dtt, My] = obj.ky_vec_dtt(dtt_type(floor(end/2) + 1));
                    [kz_vec_dtt, Mz] = obj.kz_vec_dtt(dtt_type(end));

                    % define the wavenumber based on the wavenumber
                    % components
                    k = zeros(obj.Nx, obj.Ny, obj.Nz);
                    k = bsxfun(@plus, (reshape(kx_vec_dtt, [], 1, 1).^2), k);
                    k = bsxfun(@plus, (reshape(ky_vec_dtt, 1, [], 1).^2), k);
                    k = bsxfun(@plus, (reshape(kz_vec_dtt, 1, 1, []).^2), k);
                    k = sqrt(k);
                    
                    % define product of implied period
                    M = Mx*My*Mz;
                    
            end            

        end
        
        % compute the DTT wavenumber vector in the x-direction
        function [kx_vec_dtt, M] = kx_vec_dtt(obj, dtt_type)     
            [kx_vec_dtt, M] = obj.makeDTTDim(obj.Nx, obj.dx, dtt_type);
        end
            
        % compute the DTT wavenumber vector in the x-direction
        function [ky_vec_dtt, M] = ky_vec_dtt(obj, dtt_type)     
            [ky_vec_dtt, M] = obj.makeDTTDim(obj.Ny, obj.dy, dtt_type);
        end
        
        % compute the DTT wavenumber vector in the x-direction
        function [kz_vec_dtt, M] = kz_vec_dtt(obj, dtt_type)     
            [kz_vec_dtt, M] = obj.makeDTTDim(obj.Nz, obj.dz, dtt_type);
        end
        
        % calculate highest prime factors
        function highest_prime_factors = highest_prime_factors(obj, axisymmetric)
            if nargin == 2
                switch axisymmetric
                    case 'WSWA'
                        highest_prime_factors = [max(factor(obj.Nx)), max(factor(obj.Ny * 4)), max(factor(obj.Nz))];
                    case 'WSWS'
                        highest_prime_factors = [max(factor(obj.Nx)), max(factor(obj.Ny * 2 - 2)), max(factor(obj.Nz))];
                    otherwise
                        error('Unknown axisymmetric symmetry.');
                end
            else
                highest_prime_factors = [max(factor(obj.Nx)), max(factor(obj.Ny)), max(factor(obj.Nz))];
            end
        end        
 
    end
    
    % functions that can only be accessed by class members
    methods (Access = 'protected', Static = true) 
        
        % subfunction to create the grid parameters for a single spatial
        % direction 
        function kx_vec = makeDim(Nx, dx)

            % define the discretisation of the spatial dimension such that
            % there is always a DC component
            if rem(Nx, 2) == 0
                % grid dimension has an even number of points
                nx = ((-Nx/2:Nx/2-1)/Nx).';
            else
                % grid dimension has an odd number of points
                nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx).';
            end

            % force middle value to be zero in case 1/Nx is a recurring
            % number and the series doesn't give exactly zero
            nx(floor(Nx/2) + 1) = 0;
            
            % define the wavenumber vector components
            kx_vec = (2*pi/dx) .* nx;       
        end
        
        % subfunction to create the DTT grid parameters for a single
        % spatial direction 
        function [kx_vec, M] = makeDTTDim(Nx, dx, dtt_type)
        
            % compute the implied period of the input function
            switch dtt_type
                case 1
                    M = 2 .* (Nx - 1);
                case 5
                    M = 2 .* (Nx + 1);
                otherwise
                    M = 2 .* Nx;
            end

            % calculate the wavenumbers
            switch dtt_type
                case 1
                    
                    % whole-wavenumber DTT 
                    % WSWS / DCT-I
                    n = (0:1:M/2).';
                    kx_vec = 2 .* pi .* n ./ (M .* dx);
                    
                case 2
                    
                    % whole-wavenumber DTT
                    % HSHS / DCT-II
                    n = (0:1:(M/2 - 1)).';
                    kx_vec = 2 .* pi .* n ./ (M .* dx);
                    
                case 5
                    
                    % whole-wavenumber DTT
                    % WAWA / DST-I
                    n = (1:1:(M/2 - 1)).';
                    kx_vec = 2 .* pi .* n ./ (M .* dx);
                    
                case 6
                    
                    % whole-wavenumber DTT
                    % HAHA / DST-II
                    n = (1:1:M/2).';
                    kx_vec = 2 .* pi .* n ./ (M .* dx);
                    
                case {3, 4, 7, 8}
                    
                    % half-wavenumber DTTs
                    % WSWA / DCT-III
                    % HSHA / DCT-IV
                    % WAWS / DST-III
                    % HAHS / DST-IV
                    n = (0:1:(M/2 - 1)).';
                    kx_vec = 2 .* pi .* (n + 0.5) ./ (M .* dx);
                    
            end
            
        end
        
    end
    
    % functions for non-uniform grids
    methods
        
        % subfunction to return plaid xn matrix
        function xn = get.xn(obj)
            if obj.nonuniform
                % duplicate xn vector
                switch obj.dim
                    case 1
                        xn = obj.xn_vec;
                    case 2
                        xn = repmat(obj.xn_vec, [1, obj.Ny]);
                    case 3
                        xn = repmat(obj.xn_vec, [1, obj.Ny, obj.Nz]);
                end
            else
                xn = 0;
            end
        end
        
        % subfunction to return plaid yn matrix
        function yn = get.yn(obj)
            if obj.nonuniform            
                % rotate yn vector to y direction then duplicate
                switch obj.dim
                    case 1
                        yn = 0;
                    case 2
                        yn = repmat(obj.yn_vec.', [obj.Nx, 1]);
                    case 3
                        yn = repmat(obj.yn_vec.', [obj.Nx, 1, obj.Nz]);
                end
            else
                yn = 0;
            end
        end     
        
        % subfunction to return plaid zn matrix
        function zn = get.zn(obj)
            if obj.nonuniform             
                % permute kz vector to z direction then duplicate
                switch obj.dim
                    case 1
                        zn = 0;
                    case 2
                        zn = 0;
                    case 3
                        zn = repmat(permute(obj.zn_vec, [2 3 1]), [obj.Nx, obj.Ny, 1]);
                end
            else
                zn = 0;
            end
        end   
        
        % subfunction to set non-uniform grid parameters in specified
        % dimension
        function obj = setNUGrid(obj, dim, xn_vec, dxudxn, xn_vec_sgx, dxudxn_sgx)
            
            % check the dimension to set the nonuniform grid is appropriate
            if dim > obj.dim
                error(['Cannot set nonuniform parameters for dimension ' num2str(dim) ' of ' num2str(obj.dim) '-dimensional grid.']);
            end
            
            % force non-uniform grid spacing to be column vectors, and the
            % gradients to be in the correct direction for use with bsxfun
            switch dim
                case 1
                    obj.xn_vec     = reshape(xn_vec,     [], 1);
                    obj.xn_vec_sgx = reshape(xn_vec_sgx, [], 1);
                    
                    obj.dxudxn     = reshape(dxudxn,     [], 1);
                    obj.dxudxn_sgx = reshape(dxudxn_sgx, [], 1);
                case 2
                    obj.yn_vec     = reshape(xn_vec,     [], 1);
                    obj.yn_vec_sgy = reshape(xn_vec_sgx, [], 1);
                    
                    obj.dyudyn     = reshape(dxudxn,     1, []);
                    obj.dyudyn_sgy = reshape(dxudxn_sgx, 1, []);
                case 3
                    obj.zn_vec     = reshape(xn_vec,     [], 1);
                    obj.zn_vec_sgz = reshape(xn_vec_sgx, [], 1);
                    
                    obj.dzudzn     = reshape(dxudxn,     1, 1, []);
                    obj.dzudzn_sgz = reshape(dxudxn_sgx, 1, 1, []);
            end
            
            % set non-uniform flag
            obj.nonuniform = true;
            
        end
    end
end