function sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDERAS Axisymmetric time-domain simulation of wave propagation.
%
% DESCRIPTION:
%     kspaceFirstOrderAS simulates the time-domain propagation of
%     compressional waves through an axisymmetric homogeneous or
%     heterogeneous acoustic medium. The code is functionally very similar
%     to kspaceFirstOrder2D. However, a 2D axisymmetric coordinate system
%     is used instead of a 2D Cartesian coordinate system. In this case, x
%     corresponds to the axial dimension, and y corresponds to the radial
%     dimension. In the radial dimension, the first grid point corresponds
%     to the grid origin, i.e., y = 0. In comparison, for
%     kspaceFirstOrder2D, the Cartesian point y = 0 is in the middle of the
%     computational grid.
%
%     The input structures kgrid, medium, source, and sensor are defined in
%     exactly the same way as for kspaceFirstOrder2D. However,
%     computationally, there are several key differences. First, the
%     axisymmetric code solves the coupled first-order equations accounting
%     for viscous absorption (not power law), so only medium.alpha_power =
%     2 is supported. This value is set by default, and doesn't need to be
%     defined. This also means that medium.alpha_mode and
%     medium.alpha_filter are not supported. Second, for a homogeneous
%     medium, the k-space correction used to counteract the numerical
%     dispersion introduced by the finite-difference time step is not exact
%     (as it is for the other fluid codes). However, the approximate
%     k-space correction still works very effectively, so dispersion errors
%     should still be small. See kspaceFirstOrder2D for additional details
%     on the function inputs.
%
%     In the x-dimension (axial), the FFT is used to compute spatial
%     gradients. In the y-dimension (radial), two choices of symmetry are
%     possible. These are whole-sample-symmetric on the interior radial
%     boundary (y = 0) and either whole-sample-symmetric or
%     whole-sample-asymmetric on the exterior radial boundary. These are
%     abbreviated WSWA and WSWS. The WSWA and WSWS symmetries are
%     implemented using both discrete trigonometric transforms (DTTs), and
%     via the FFT by manually mirroring the domain. The latter options are
%     abbreviated as WSWA-FFT and WSWS-FFT. The WSWA/WSWS options and the
%     corresponding WSWA-FFT/WSWS-FFT options agree to machine precision.
%     When using the PML, the choice of symmetry doesn't matter, and all
%     options give very similar results (to several decimal places).
%     Computationally, the DTT implementations are more efficient, but
%     require additional compiled MATLAB functions (not currently part of
%     k-Wave). The symmetry can be set by using the optional input
%     'RadialSymmetry'. The WSWA-FFT symmetry is set by default. 
%        
% USAGE:
%     sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor)
%     sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
% The minimum fields that must be assigned to run an initial value problem
% (for example, a photoacoustic forward simulation) are marked with a *. 
%
%     kgrid*                 - k-Wave grid object returned by kWaveGrid
%                              containing axisymmetric and k-space grid
%                              fields
%     kgrid.t_array*         - evenly spaced array of time values [s] (set
%                              to 'auto' by kWaveGrid) 
%
%
%     medium.sound_speed*    - sound speed distribution within the acoustic
%                              medium [m/s] 
%     medium.sound_speed_ref - reference sound speed used within the
%                              k-space operator (phase correction term)
%                              [m/s]
%     medium.density*        - density distribution within the acoustic
%                              medium [kg/m^3] 
%     medium.BonA            - parameter of nonlinearity
%     medium.alpha_coeff     - power law absorption coefficient 
%                              [dB/(MHz^y cm)] 
%     medium.alpha_sign      - scalar used to control the sign of the
%                              absorption term in the equation of state 
%
%     source.p0*             - initial pressure within the acoustic medium
%     source.p               - time varying pressure at each of the source
%                              positions given by source.p_mask 
%     source.p_mask          - binary matrix specifying the positions of
%                              the time varying pressure source
%                              distribution 
%     source.p_mode          - optional input to control whether the input
%                              pressure is injected as a mass source or
%                              enforced as a dirichlet boundary condition;
%                              valid inputs are 'additive' (the default) or
%                              'dirichlet'    
%     source.ux              - time varying particle velocity in the
%                              x-direction (axial) at each of the source
%                              positions given by source.u_mask
%     source.uy              - time varying particle velocity in the
%                              y-direction (radial) at each of the source
%                              positions given by source.u_mask  
%     source.u_mask          - binary matrix specifying the positions of
%                              the time varying particle velocity
%                              distribution
%     source.u_mode          - optional input to control whether the input
%                              velocity is applied as a force source or
%                              enforced as a dirichlet boundary condition;
%                              valid inputs are 'additive' (the default) or
%                              'dirichlet'
%
%     sensor.mask*           - binary matrix or a set of Cartesian points
%                              where the pressure is recorded at each
%                              time-step  
%     sensor.record          - cell array of the acoustic parameters to
%                              record in the form sensor.record = {'p',
%                                'u', ...}; valid inputs are:  
%                                'p' (acoustic pressure)
%                                'p_max' (maximum pressure)
%                                'p_min' (minimum pressure)
%                                'p_rms' (RMS pressure)
%                                'p_final' (final pressure field at all grid points)
%                                'p_max_all' (maximum pressure at all grid points)
%                                'p_min_all' (minimum pressure at all grid points)
%                                'u' (particle velocity)
%                                'u_max' (maximum particle velocity)
%                                'u_min' (minimum particle velocity)
%                                'u_rms' (RMS particle velocity)
%                                'u_final' (final particle velocity field at all grid points)
%                                'u_max_all' (maximum particle velocity at all grid points)
%                                'u_min_all' (minimum particle velocity at all grid points)
%                                'u_non_staggered' (particle velocity on non-staggered grid)
%                                'I' (time varying acoustic intensity)
%                                'I_avg' (average acoustic intensity) 
%     sensor.record_start_index 
%                            - time index at which the sensor should start
%                              recording the data specified by 
%                              sensor.record (default = 1) 
%     sensor.time_reversal_boundary_data 
%                            - time varying pressure enforced as a
%                              Dirichlet boundary condition over
%                              sensor.mask
%     sensor.frequency_response 
%                            - two element array specifying the center
%                              frequency and percentage bandwidth of a
%                              frequency domain Gaussian filter applied to
%                              the sensor_data
%
% Note: For heterogeneous medium parameters, medium.sound_speed and
% medium.density must be given in matrix form with the same dimensions as
% kgrid. For homogeneous medium parameters, these can be given as single
% numeric values. If the medium is homogeneous and velocity inputs or
% outputs are not required, it is not necessary to specify medium.density.
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings.
%
%     'CartInterp'           - Interpolation mode used to extract the
%                              pressure when a Cartesian sensor mask is
%                              given. If set to 'nearest' and more than one
%                              Cartesian point maps to the same grid point,
%                              duplicated data points are discarded and
%                              sensor_data will be returned with less
%                              points than that specified by sensor.mask
%                              (default = 'linear'). 
%     'CreateLog'            - Boolean controlling whether the command line
%                              output is saved using the diary function
%                              with a date and time stamped filename
%                              (default = false).
%     'DataCast'             - String input of the data type that variables
%                              are cast to before computation. For example,
%                              setting to 'single' will speed up the
%                              computation time (due to the improved
%                              efficiency of fftn and ifftn for this data
%                              type) at the expense of a loss in precision. 
%                              This variable is also useful for utilising
%                              GPU parallelisation through libraries such
%                              as the Parallel Computing Toolbox by setting
%                              'DataCast' to 'gpuArray-single' (default =
%                              'off').
%     'DataRecast'           - Boolean controlling whether the output data
%                              is cast back to double precision. If set to
%                              false, sensor_data will be returned in the
%                              data format set using the 'DataCast' option.
%     'DisplayMask'          - Binary matrix overlaid onto the animated
%                              simulation display. Elements set to 1 within
%                              the display mask are set to black within the
%                              display (default = sensor.mask).
%     'HDFCompressionLevel'  - Compression level used for writing the input
%                              HDF5 file when using 'SaveToDisk' or
%                              kspaceFirstOrderASC. Can be set to an
%                              integer between 0 (no compression, the
%                              default) and 9 (maximum compression). The
%                              compression is lossless. Increasing the
%                              compression level will reduce the file size
%                              if there are portions of the medium that are
%                              homogeneous, but will also increase the time
%                              to create the HDF5 file.
%     'LogScale'             - Boolean controlling whether the pressure
%                              field is log compressed before display
%                              (default = false). The data is compressed by
%                              scaling both the positive and negative
%                              values between 0 and 1 (truncating the data
%                              to the given plot scale), adding a scalar
%                              value (compression factor) and then using
%                              the corresponding portion of a log10 plot
%                              for the compression (the negative parts are
%                              remapped to be negative thus the default
%                              color scale will appear unchanged). The
%                              amount of compression can be controlled by
%                              adjusting the compression factor which can
%                              be given in place of the Boolean input. The
%                              closer the compression factor is to zero,
%                              the steeper the corresponding part of the 
%                              log10 plot used, and the greater the
%                              compression (the default compression factor
%                              is 0.02).
%     'MeshPlot'             - Boolean controlling whether mesh is used in 
%                              place of imagesc to plot the pressure field
%                              (default = false). When 'MeshPlot' is set to
%                              true, the default display mask is set to
%                              'off'.  
%     'MovieArgs'            - Settings for VideoWriter. Parameters must be
%                              given as {'param', value, ...} pairs within
%                              a cell array (default = {}), where 'param'
%                              corresponds to a writable property of a
%                              VideoWriter object. 
%     'MovieName'            - Name of the movie produced when
%                              'RecordMovie' is set to true (default =
%                              'date-time-kspaceFirstOrder2D'). 
%     'MovieProfile'         - Profile input passed to VideoWriter.
%     'PlotFreq'             - The number of iterations which must pass 
%                              before the simulation plot is updated
%                              (default = 10).
%     'PlotLayout'           - Boolean controlling whether a four panel
%                              plot of the initial simulation layout is
%                              produced (initial pressure, sensor mask,
%                              sound speed, density) (default = false).
%     'PlotPML'              - Boolean controlling whether the perfectly
%                              matched layer is shown in the simulation
%                              plots. If set to false, the PML is not
%                              displayed (default = true). 
%     'PlotScale'            - [min, max] values used to control the
%                              scaling for imagesc (visualisation). If set
%                              to 'auto', a symmetric plot scale is chosen
%                              automatically for each plot frame.
%     'PlotSim'              - Boolean controlling whether the simulation
%                              iterations are progressively plotted
%                              (default = true). 
%     'PMLAlpha'             - Absorption within the perfectly matched 
%                              layer in Nepers per grid point (default =
%                              2).
%     'PMLInside'            - Boolean controlling whether the perfectly 
%                              matched layer is inside or outside the grid.
%                              If set to false, the input grids are
%                              enlarged by PMLSize before running the
%                              simulation (default = true).
%     'PMLSize'              - Size of the perfectly matched layer in grid
%                              points. By default, the PML is added evenly
%                              to all sides of the grid, however, both
%                              PMLSize and PMLAlpha can be given as two
%                              element arrays to specify the x and y
%                              properties, respectively. To remove the PML,
%                              set the appropriate PMLAlpha to zero rather
%                              than forcing the PML to be of zero size
%                              (default = 20). 
%     'RadialSymmetry'       - Radial symmetry assumed at the inner and
%                              outer domain boundaries by the simulation,
%                              where W: whole sample, H: half sample, S:
%                              symmetric, A: antisymmetric. Valid inputs
%                              are 'WSWA-FFT', 'WSWS-FFT', 'WSWA', and
%                              'WSWS'. The first two options are computed
%                              using the FFT by mirroring the domain
%                              appropriately. The final two options are
%                              computed using the discrete cosine and sine
%                              transformation and require additional
%                              functions (not currently part of k-Wave).
%     'RecordMovie'          - Boolean controlling whether the displayed
%                              image frames are captured and stored as a
%                              movie using VideoWriter (default = false).
%     'Smooth'               - Boolean controlling whether source.p0,
%                              medium.sound_speed, and medium.density are
%                              smoothed using smooth before computation.
%                              'Smooth' can either be given as a single
%                              Boolean value or as a 3 element array to
%                              control the smoothing of source.p0,
%                              medium.sound_speed, and medium.density,
%                              independently (default = [true, false,
%                              false]). 
%
% OUTPUTS:
% If sensor.record is not defined by the user:
%     sensor_data            - time varying pressure recorded at the sensor
%                              positions given by sensor.mask
%
% If sensor.record is defined by the user:
%     sensor_data.p          - time varying pressure recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'p' is set)
%     sensor_data.p_max      - maximum pressure recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'p_max' is set)  
%     sensor_data.p_min      - minimum pressure recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'p_min' is set)  
%     sensor_data.p_rms      - rms of the time varying pressure recorded at
%                              the sensor positions given by sensor.mask
%                              (returned if 'p_rms' is set)
%     sensor_data.p_final    - final pressure field at all grid points
%                              within the domain (returned if 'p_final' is
%                              set)
%     sensor_data.p_max_all  - maximum pressure recorded at all grid points
%                              within the domain (returned if 'p_max_all'
%                              is set) 
%     sensor_data.p_min_all  - minimum pressure recorded at all grid points
%                              within the domain (returned if 'p_min_all'
%                              is set)
%     sensor_data.ux         - time varying particle velocity in the
%                              x-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'u' is
%                              set) 
%     sensor_data.uy         - time varying particle velocity in the
%                              y-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'u' is
%                              set)
%     sensor_data.ux_max     - maximum particle velocity in the x-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_max' is set)     
%     sensor_data.uy_max     - maximum particle velocity in the y-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_max' is set) 
%     sensor_data.ux_min     - minimum particle velocity in the x-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_min' is set)
%     sensor_data.uy_min     - minimum particle velocity in the y-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_min' is set)   
%     sensor_data.ux_rms     - rms of the time varying particle velocity in
%                              the x-direction recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'u_rms' is set)     
%     sensor_data.uy_rms     - rms of the time varying particle velocity in
%                              the y-direction recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'u_rms' is set)        
%     sensor_data.ux_final   - final particle velocity field in the
%                              x-direction at all grid points within the
%                              domain (returned if 'u_final' is set) 
%     sensor_data.uy_final   - final particle velocity field in the
%                              y-direction at all grid points within the
%                              domain (returned if 'u_final' is set)   
%     sensor_data.ux_max_all - maximum particle velocity in the x-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_max_all' is set) 
%     sensor_data.uy_max_all - maximum particle velocity in the y-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_max_all' is set) 
%     sensor_data.ux_min_all - minimum particle velocity in the x-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_min_all' is set)
%     sensor_data.uy_min_all - minimum particle velocity in the y-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_min_all' is set) 
%     sensor_data.ux_non_staggered 
%                            - time varying particle velocity in the
%                              x-direction recorded at the sensor positions
%                              given by sensor.mask after shifting to the
%                              non-staggered grid (returned if
%                              'u_non_staggered' is set)
%     sensor_data.uy_non_staggered 
%                            - time varying particle velocity in the
%                              y-direction recorded at the sensor positions
%                              given by sensor.mask after shifting to the
%                              non-staggered grid (returned if
%                              'u_non_staggered' is set)
%     sensor_data.Ix         - time varying acoustic intensity in the
%                              x-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'I' is
%                              set)
%     sensor_data.Iy         - time varying acoustic intensity in the
%                              y-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'I' is
%                              set)
%     sensor_data.Ix_avg     - average acoustic intensity in the
%                              x-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'I_avg' is
%                              set)  
%     sensor_data.Iy_avg     - average acoustic intensity in the
%                              y-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'I_avg' is
%                              set)  
%
% ABOUT:
%     author                 - Bradley Treeby, Elliott Wise, and Ben Cox
%     date                   - 28th September 2017
%     last update            - 24th March 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Bradley Treeby, Elliott Wise, and Ben Cox
%
% See also kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kWaveGrid

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

% suppress mlint warnings that arise from using subscripts
%#ok<*NASGU>
%#ok<*COLND>
%#ok<*NODEF>
%#ok<*INUSL>

% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================

% start the timer and store the start time
start_time = clock;
tic;

% set the name of the simulation code
MFILE = mfilename;

% get the number of inputs and outputs (nargin and nargout can't be used in
% subscripts in MATLAB 2016b or later)
num_inputs  = nargin;
num_outputs = nargout;

% run subscript to check inputs
kspaceFirstOrder_inputChecking;

% =========================================================================
% CALCULATE MEDIUM PROPERTIES ON STAGGERED GRID
% =========================================================================

% interpolate the values of the density at the staggered grid locations
% where sgx = (x + dx/2, y) and sgy = (x, y + dy/2)
if numDim(rho0) == 2 && flags.use_sg
    
    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interpn(kgrid.x, kgrid.y, rho0, kgrid.x + kgrid.dx/2, kgrid.y, '*linear');
    rho0_sgy = interpn(kgrid.x, kgrid.y, rho0, kgrid.x, kgrid.y + kgrid.dy/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx));
    rho0_sgy(isnan(rho0_sgy)) = rho0(isnan(rho0_sgy));
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
    rho0_sgy = rho0;
end

% invert rho0 so it doesn't have to be done each time step
rho0_sgx_inv = 1./rho0_sgx;
rho0_sgy_inv = 1./rho0_sgy;

% clear unused variables if not using them in _saveToDisk
if ~flags.save_to_disk
    clear rho0_sgx rho0_sgy 
end

% =========================================================================
% PREPARE DERIVATIVE AND PML OPERATORS
% =========================================================================

% get the PML operators based on the reference sound speed and PML settings
pml_x     = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, false,                1, false);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, true && flags.use_sg, 1, false);
pml_y     = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, false,                2, true);
pml_y_sgy = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, true && flags.use_sg, 2, true);

% define the k-space, derivative, and shift operators
% for the x (axial) direction, the operators are the same as normal
ddx_k_shift_pos = ifftshift( 1i * kgrid.kx_vec .* exp( 1i * kgrid.kx_vec * kgrid.dx/2) );
ddx_k_shift_neg = ifftshift( 1i * kgrid.kx_vec .* exp(-1i * kgrid.kx_vec * kgrid.dx/2) );

% for the y (radial) direction
% when using DTTs:
%    - there is no explicit grid shift (this is done by choosing DTTs
%      with the appropriate symmetry)
%    - ifftshift isn't required as the wavenumbers start from DC
% when using FFTs:
%    - the grid is expanded, and the fields replicated in the radial
%      dimension to give the required symmetry
%    - the derivative and shift operators are defined as normal
switch radial_symmetry
    case {'WSWA-FFT', 'WSWS-FFT'}

        % create a new kWave grid object with expanded radial grid
        switch radial_symmetry
            case 'WSWA-FFT'

                % extend grid by a factor of x4 to account for
                % symmetries in WSWA
                kgrid_exp = kWaveGrid(kgrid.Nx, kgrid.dx, kgrid.Ny * 4, kgrid.dy);                    

            case 'WSWS-FFT'

                % extend grid by a factor of x2 - 2 to account for
                % symmetries in WSWS
                kgrid_exp = kWaveGrid(kgrid.Nx, kgrid.dx, kgrid.Ny * 2 - 2, kgrid.dy);                          

        end

        % define operators, rotating y-direction for use with bsxfun
        ddy_k       = ifftshift( 1i * kgrid_exp.ky_vec ).';
        y_shift_pos = ifftshift( exp( 1i * kgrid_exp.ky_vec * kgrid_exp.dy/2) ).';
        y_shift_neg = ifftshift( exp(-1i * kgrid_exp.ky_vec * kgrid_exp.dy/2) ).';

        % define the k-space operator
        if flags.use_kspace
            kappa = ifftshift(sinc(c_ref .* kgrid_exp.k .* kgrid.dt / 2));
            if (flags.source_p && strcmp(source.p_mode, 'additive')) || ((flags.source_ux || flags.source_uy) && strcmp(source.u_mode, 'additive'))
                source_kappa = ifftshift(cos (c_ref .* kgrid_exp.k .* kgrid.dt / 2));
            end
        else
            kappa = 1;
            source_kappa = 1;
        end

    case {'WSWA', 'WSWS'}
        
        switch radial_symmetry
            
            case 'WSWA'

                % get the wavenumbers and implied length for the DTTs
                [ky_vec, M] = kgrid.ky_vec_dtt(DCT3);

                % define the derivative operators
                ddy_k_wswa = -ky_vec.';
                ddy_k_hahs =  ky_vec.';
                
            case 'WSWS'
                
                % get the wavenumbers and implied length for the DTTs
                [ky_vec, M] = kgrid.ky_vec_dtt(DCT1);

                % define the derivative operators
                ddy_k_wsws = -ky_vec(2:end).';
                ddy_k_haha =  ky_vec(2:end).';                
                
        end

        % define the k-space operator
        if flags.use_kspace
            
            % define scalar wavenumber
            k_dtt = sqrt( repmat(ifftshift(kgrid.kx_vec).^2, 1, kgrid.Ny) + repmat((ky_vec.').^2, kgrid.Nx, 1) );
            
            % define k-space operators
            kappa = sinc(c_ref .* k_dtt .* kgrid.dt / 2);
            if (flags.source_p && strcmp(source.p_mode, 'additive')) || ((flags.source_ux || flags.source_uy) && strcmp(source.u_mode, 'additive'))
                source_kappa = cos (c_ref .* k_dtt .* kgrid.dt / 2);
            end
            
            % cleanup unused variables
            clear k_dtt;
            
        else
            kappa = 1;
            source_kappa = 1;
        end

end

% define staggered and non-staggered grid axial distance
y_vec    = (kgrid.y_vec - kgrid.y_vec(1)).';
y_vec_sg = (kgrid.y_vec - kgrid.y_vec(1) + kgrid.dy/2).';
    
% option to run simulations without the spatial staggered grid is not
% supported for the axisymmetric code 
if ~flags.use_sg
    error('Optional input ''UseSG'' is not supported for axisymmetric simulations.');
end
   
% =========================================================================
% SAVE DATA TO DISK FOR RUNNING SIMULATION EXTERNAL TO MATLAB
% =========================================================================

% save to disk option for saving the input matrices to disk for running
% simulations using k-Wave++
if flags.save_to_disk
    
    % run subscript to save files to disk
    kspaceFirstOrder_saveToDisk;
    
    % exit matlab computation if required
    if flags.save_to_disk_exit
        sensor_data = [];
        return
    end
    
end

% =========================================================================
% DATA CASTING
% =========================================================================

% preallocate the loop variables using the castZeros anonymous function
% (this creates a matrix of zeros in the data type specified by data_cast)
p                   = castZeros([kgrid.Nx, kgrid.Ny]);
p_k                 = castZeros([kgrid.Nx, kgrid.Ny]);
dpdx_sgx            = castZeros([kgrid.Nx, kgrid.Ny]);
dpdy_sgy            = castZeros([kgrid.Nx, kgrid.Ny]);
rhox                = castZeros([kgrid.Nx, kgrid.Ny]);
rhoy                = castZeros([kgrid.Nx, kgrid.Ny]);
ux_sgx              = castZeros([kgrid.Nx, kgrid.Ny]);
uy_sgy              = castZeros([kgrid.Nx, kgrid.Ny]);
duxdx               = castZeros([kgrid.Nx, kgrid.Ny]);
duydy               = castZeros([kgrid.Nx, kgrid.Ny]);

% preallocate expanded grid variables needed if using the FFT versions
if any(strcmp(radial_symmetry, {'WSWA-FFT', 'WSWS-FFT'}))
    p_exp           = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    ux_exp          = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    uy_exp          = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    uy_on_y         = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    dpdx_sgx_exp    = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    dpdy_sgy_exp    = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    duxdx_exp       = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);
    duydy_exp       = castZeros([kgrid_exp.Nx, kgrid_exp.Ny]);    
end

% run subscript to cast the remaining loop variables to the data type
% specified by data_cast 
if ~strcmp(data_cast, 'off')
    kspaceFirstOrder_dataCast;
end

% =========================================================================
% CREATE INDEX VARIABLES
% =========================================================================

% setup the time index variable
if ~flags.time_rev
    index_start = 1;
    index_step  = 1;
    index_end   = kgrid.Nt; 
else
    
    % reverse the order of the input data
    sensor.time_reversal_boundary_data = fliplr(sensor.time_reversal_boundary_data);
    index_start = 1;
    index_step = 1;
    
    % stop one time point before the end so the last points are not
    % propagated
    index_end = kgrid.Nt - 1;    
    
end

% =========================================================================
% PREPARE VISUALISATIONS
% =========================================================================

% pre-compute suitable axes scaling factor
if flags.plot_layout || flags.plot_sim
    [x_sc, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec]));  %#ok<ASGLU>
end

% run subscript to plot the simulation layout if 'PlotLayout' is set to true
if flags.plot_layout
    kspaceFirstOrder_plotLayout;
end

% initialise the figure used for animation if 'PlotSim' is set to 'true'
if flags.plot_sim
    kspaceFirstOrder_initialiseFigureWindow;
end 

% initialise movie parameters if 'RecordMovie' is set to 'true'
if flags.record_movie
    kspaceFirstOrder_initialiseMovieParameters;
end

% =========================================================================
% LOOP THROUGH TIME STEPS
% =========================================================================

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
disp('  starting time loop...');

% restart timing variables
loop_start_time = clock;
tic;

% start time loop
for t_index = index_start:index_step:index_end

    % enforce time reversal bounday condition
    if flags.time_rev
      
        % load pressure value and enforce as a Dirichlet boundary condition
        p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, t_index);

        % compute rhox and rhoy using an adiabatic equation of state
        rhox_mod = 0.5 .* p ./ (c0.^2);
        rhoy_mod = 0.5 .* p ./ (c0.^2);
        rhox(sensor_mask_index) = rhox_mod(sensor_mask_index);
        rhoy(sensor_mask_index) = rhoy_mod(sensor_mask_index);
            
    end
 
    % calculate dp/dx and dp/dy at the current time step
    switch radial_symmetry          
        case {'WSWA-FFT', 'WSWS-FFT'}
            
            % mirror p in the radial dimension using appropriate symmetry
            switch radial_symmetry
                case 'WSWA-FFT'
                    p_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         p;
                    p_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2)     = -fliplr(p(:, 2:end));
                    p_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)     =        -p;
                    p_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4)     =  fliplr(p(:, 2:end));
                case 'WSWS-FFT'
                    p_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         p;
                    p_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2 - 2) =  fliplr(p(:, 2:end - 1));
            end
            
            % compute gradients using k-space pseudospectral scheme
            p_k = kappa .* fft2(p_exp);
            dpdx_sgx_exp = real(ifft2( bsxfun(@times, ddx_k_shift_pos, p_k) ));
            dpdy_sgy_exp = real(ifft2( bsxfun(@times, ddy_k .* y_shift_pos, p_k) ));

            % trim the gradients
            dpdx_sgx = dpdx_sgx_exp(:, 1:kgrid.Ny);
            dpdy_sgy = dpdy_sgy_exp(:, 1:kgrid.Ny);

        case 'WSWA'

            % compute forward transform in the radial direction using
            % DCT3 (WSWA symmetry), with transforms in the axial
            % direction computed using the FFT 
            p_k = kappa .* fft(dtt1D(p, DCT3, 2), [], 1);

            % for the axial derivative, compute inverse transform in
            % the radial direction using DCT3^-1 = DCT2 (symmetry and
            % grid staggering don't change)
            dpdx_sgx = dtt1D(real(ifft(bsxfun(@times, ddx_k_shift_pos, p_k), [], 1)), DCT2, 2) ./ M;

            % for the radial derivative, compute inverse transform in
            % the radial direction using DST4^-1 = DST4 (symmetry
            % changes to HAHS as output is differentiated and moved to
            % staggered grid)  
            dpdy_sgy = dtt1D(real(ifft(bsxfun(@times, ddy_k_wswa, p_k), [], 1)), DST4, 2) ./ M;

        case 'WSWS'

            % compute forward transform in the radial direction using
            % DCT1 (WSWS symmetry), with transforms in the axial
            % direction computed using the FFT 
            p_k = kappa .* fft(dtt1D(p, DCT1, 2), [], 1);

            % for the axial derivative, compute inverse transform in
            % the radial direction using DCT1^-1 = DCT1 (symmetry and
            % grid staggering don't change)
            dpdx_sgx = dtt1D(real(ifft(bsxfun(@times, ddx_k_shift_pos, p_k), [], 1)), DCT1, 2) ./ M;

            % for the radial deriviative, remove the left end point,
            % then compute inverse transform using DST2^-1 = DST3
            % (symmetry changes to HAHA as output is differentiated and
            % moved to staggered grid)
            dpdy_sgy(:, 1:end - 1) = dtt1D(real(ifft(bsxfun(@times, ddy_k_wsws, p_k(:, 2:end)), [], 1)), DST3, 2) ./ M;

            % mirror the final point according to HAHA symmetry to
            % maintain the size of the matrix
            dpdy_sgy(:, end) = -dpdy_sgy(:, end - 1);

    end
    
    % calculate ux and uy at the next time step using dp/dx and dp/dy at
    % the current time step
    ux_sgx = bsxfun(@times, pml_x_sgx, ...
        bsxfun(@times, pml_x_sgx, ux_sgx) ...
        - dt .* rho0_sgx_inv .* dpdx_sgx ...
        );
    uy_sgy = bsxfun(@times, pml_y_sgy, ...
        bsxfun(@times, pml_y_sgy, uy_sgy) ...
        - dt .* rho0_sgy_inv .* dpdy_sgy ...
        );  
    
    % add in the velocity source terms
    if flags.source_ux >= t_index
        switch source.u_mode
            case 'dirichlet'
            
                % enforce the source values as a dirichlet boundary condition
                ux_sgx(u_source_pos_index) = source.ux(u_source_sig_index, t_index);
            
            case 'additive'
                
                % extract the source values into a matrix
                source_mat = castZeros(size(ux_sgx));
                source_mat(u_source_pos_index) = source.ux(u_source_sig_index, t_index);

                % apply the k-space correction
                switch radial_symmetry
                    case {'WSWA-FFT', 'WSWS-FFT'}
                        
                        % mirror the source matrix in the radial dimension
                        % using the appropriate symmetry 
                        switch radial_symmetry
                            case 'WSWA-FFT'
                                
                                % mirror ux source matrix using WSWA symmetry 
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)  =         source_mat;
                                source_mat_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2)  = -fliplr(source_mat(:, 2:end));
                                source_mat_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)  =        -source_mat;
                                source_mat_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4)  =  fliplr(source_mat(:, 2:end));
                                
                            case 'WSWS-FFT'
                                
                                % mirror ux source matrix using WSWS symmetry 
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)      =         source_mat;
                                source_mat_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2 - 2)  =  fliplr(source_mat(:, 2:end - 1));
                                
                        end
                        
                        % apply source correction
                        source_mat_exp = real(ifft2(source_kappa .* fft2(source_mat_exp)));
                        
                        % trim source to correct size
                        source_mat = source_mat_exp(:, 1:kgrid.Ny);
                        
                    case 'WSWA'
                        
                        % compute forward transform in the radial direction
                        % using DCT3 and inverse transform using DCT3^-1 =
                        % DCT2 (WSWA symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DCT3, 2), [], 1)...
                            , [], 1)), DCT2, 2) ./ M;                                         
                        
                    case 'WSWS'
                        
                        % compute forward transform in the radial direction
                        % using DCT1 and inverse transform using DCT1^-1 =
                        % DCT1 (WSWS symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DCT1, 2), [], 1)...
                            , [], 1)), DCT1, 2) ./ M;                            
                        
                end

                % add the source values to the existing field values,
                % including the k-space correction
                ux_sgx = ux_sgx + source_mat;
                
            case 'additive-no-correction'
            
                % add the source values to the existing field values 
                ux_sgx(u_source_pos_index) = ux_sgx(u_source_pos_index) + source.ux(u_source_sig_index, t_index);
            
        end
    end
    if flags.source_uy >= t_index
        switch source.u_mode
            case 'dirichlet'
            
                % enforce the source values as a dirichlet boundary condition        
                uy_sgy(u_source_pos_index) = source.uy(u_source_sig_index, t_index);
            
            case 'additive'
                
                % extract the source values into a matrix
                source_mat = castZeros(size(uy_sgy));
                source_mat(u_source_pos_index) = source.uy(u_source_sig_index, t_index);

                % apply the k-space correction
                switch radial_symmetry
                    case {'WSWA-FFT', 'WSWS-FFT'}
                        
                        % mirror the source matrix in the radial dimension
                        % using the appropriate symmetry 
                        switch radial_symmetry
                            case 'WSWA-FFT'
                                
                                % mirror uy source matrix using HAHS symmetry 
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)  =         source_mat;
                                source_mat_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2)  =  fliplr(source_mat);                
                                source_mat_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)  =        -source_mat;
                                source_mat_exp(:, kgrid.Ny*3 + 1:kgrid.Ny*4)  = -fliplr(source_mat);   
                                
                            case 'WSWS-FFT'
                                
                                % mirror uy source matrix using HAHA symmetry
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1 - 1)  =         source_mat(:, 1:end - 1);
                                source_mat_exp(:, kgrid.Ny*1 + 0:kgrid.Ny*2 - 2)  = -fliplr(source_mat(:, 1:end - 1));

                        end
                        
                        % apply source correction
                        source_mat_exp = real(ifft2(source_kappa .* fft2(source_mat_exp)));
                        
                        % trim source to correct size
                        source_mat = source_mat_exp(:, 1:kgrid.Ny);
                        
                    case 'WSWA'
                        
                        % compute forward transform in the radial direction
                        % using DST4 and inverse transform using DST4^-1 =
                        % DST4 (HAHS symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DST4, 2), [], 1)...
                            , [], 1)), DST4, 2) ./ M;                                         
                        
                    case 'WSWS'
                        
                        % compute forward transform in the radial direction
                        % using DST2 and inverse transform using DST2^-1 =
                        % DST3 (HAHA symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DST2, 2), [], 1)...
                            , [], 1)), DST3, 2) ./ M;                            
                        
                end

                % add the source values to the existing field values,
                % including the k-space correction
                uy_sgy = uy_sgy + source_mat;
                
            case 'additive-no-correction'
            
                % add the source values to the existing field values 
                uy_sgy(u_source_pos_index) = uy_sgy(u_source_pos_index) + source.uy(u_source_sig_index, t_index);
            
        end
    end
        
    % calculate dux/dx and duy/dy at the next time step
    switch radial_symmetry
        case {'WSWA-FFT', 'WSWS-FFT'}
            
            % mirror data in radial dimension using appropriate symmetry
            switch radial_symmetry
                case 'WSWA-FFT'

                    % mirror ux_sgx using WSWA symmetry 
                    ux_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)  =         ux_sgx;
                    ux_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2)  = -fliplr(ux_sgx(:, 2:end));
                    ux_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)  =        -ux_sgx;
                    ux_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4)  =  fliplr(ux_sgx(:, 2:end));

                    % mirror uy_sgy using HAHS symmetry 
                    uy_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)  =         uy_sgy;
                    uy_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2)  =  fliplr(uy_sgy);                
                    uy_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)  =        -uy_sgy;
                    uy_exp(:, kgrid.Ny*3 + 1:kgrid.Ny*4)  = -fliplr(uy_sgy);            

                    % mirror uy/y in the radial dimension using HSHA symmetry
                    uy_on_y(:, kgrid.Ny*0 + 1:kgrid.Ny*1) =         bsxfun(@times, 1./y_vec_sg, uy_sgy);
                    uy_on_y(:, kgrid.Ny*1 + 1:kgrid.Ny*2) = -fliplr(bsxfun(@times, 1./y_vec_sg, uy_sgy));
                    uy_on_y(:, kgrid.Ny*2 + 1:kgrid.Ny*3) =        -bsxfun(@times, 1./y_vec_sg, uy_sgy);
                    uy_on_y(:, kgrid.Ny*3 + 1:kgrid.Ny*4) =  fliplr(bsxfun(@times, 1./y_vec_sg, uy_sgy)); 
                    
                case 'WSWS-FFT'

                    % mirror ux_sgx using WSWS symmetry 
                    ux_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)      =         ux_sgx;
                    ux_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2 - 2)  =  fliplr(ux_sgx(:, 2:end - 1));

                    % mirror uy_sgy using HAHA symmetry
                    uy_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1 - 1)  =         uy_sgy(:, 1:end - 1);
                    uy_exp(:, kgrid.Ny*1 + 0:kgrid.Ny*2 - 2)  = -fliplr(uy_sgy(:, 1:end - 1));

                    % mirror uy/y in the radial dimension using HSHS symmetry
                    uy_on_y(:, kgrid.Ny*0 + 1:kgrid.Ny*1 - 1) =        bsxfun(@times, 1./y_vec_sg(1:end - 1), uy_sgy(:, 1:end - 1));
                    uy_on_y(:, kgrid.Ny*1 + 0:kgrid.Ny*2 - 2) = fliplr(bsxfun(@times, 1./y_vec_sg(1:end - 1), uy_sgy(:, 1:end - 1)));
                    
            end
                    
            % compute gradients using k-space pseudospectral scheme
            duxdx_exp = real(ifft2( bsxfun(@times, ddx_k_shift_neg, kappa .* fft2(ux_exp)) ));
            duydy_exp = real(ifft2( kappa .* bsxfun(@times, y_shift_neg, ...
                bsxfun(@times, ddy_k, fft2(uy_exp)) + fft2(uy_on_y)) ...
                ));                

            % trim the gradients
            duxdx = duxdx_exp(:, 1:kgrid.Ny);
            duydy = duydy_exp(:, 1:kgrid.Ny);

        case 'WSWA'

            % for the axial derivative, compute forward transform in
            % the radial direction using DCT3 (WSWA symmetry), and the
            % inverse transform using DCT3^-1 = DST2 (symmetry and grid
            % staggering don't change) 
            duxdx = dtt1D(real(ifft(...
                        kappa .* bsxfun(@times, ddx_k_shift_neg, fft(dtt1D(ux_sgx, DCT3, 2), [], 1)) ...
                        , [], 1)), DCT2, 2) ./ M;

            % for the radial derivative, compute forward transform of
            % uy using DST4 (HAHS symmetry) and of uy/y using DCT4
            % (HSHA symmetry), and the inverse transform using DCT3^-1
            % = DCT2 (symmetry changes to WSWA as output is
            % differentiated and moved to regular grid)                    
            duydy = dtt1D(real(ifft(kappa .* (...
                        bsxfun(@times, ddy_k_hahs, fft(dtt1D(uy_sgy, DST4, 2), [], 1)) + ...
                        fft(dtt1D(bsxfun(@times, 1./y_vec_sg, uy_sgy), DCT4, 2), [], 1) ...
                        ), [], 1)), DCT2, 2) ./ M;                    

        case 'WSWS'

            % for the axial derivative, compute forward transform in
            % the radial direction using DCT1 (WSWS symmetry), and the
            % inverse transform using DCT1^-1 = DST1 (symmetry and grid
            % staggering don't change), with transforms in the axial
            % direction computed using the FFT
            duxdx = dtt1D(real(ifft(...
                        kappa .* bsxfun(@times, ddx_k_shift_neg, fft(dtt1D(ux_sgx, DCT1, 2), [], 1)) ...
                        , [], 1)), DCT1, 2) ./ M;

            % for the radial derivative, compute forward transform of
            % uy using DST2 (HAHA symmetry) adding a zero endpoint, with
            % transforms in the axial direction computed using the FFT
            duydy(:, 1) = 0;
            duydy(:, 2:end) = bsxfun(@times, ddy_k_haha, fft(dtt1D(uy_sgy(:, 1:end - 1), DST2, 2), [], 1));
            
            % compute forward transform of uy/y using DCT2 (HSHS symmetry),
            % with transforms in the axial direction computed using the
            % FFT, and add to duydy
            duydy(:, 1:end-1) = duydy(:, 1:end-1) + fft(dtt1D(bsxfun(@times, 1./y_vec_sg(:, 1:end - 1), uy_sgy(:, 1:end - 1)), DCT2, 2), [], 1);
            
            % apply k-space operator, and take inverse transform using
            % DCT1^-1 = DCT1 (symmetry changes to WSWS as output is
            % differentiated and moved to regular grid) 
            duydy = dtt1D(real(ifft(kappa .* duydy, [], 1)), DCT1, 2) ./ M;

    end
    
    % calculate rhox and rhoy at the next time step
    if ~flags.nonlinear
        
        % use linearised mass conservation equation
        rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0 .* duxdx);
        rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0 .* duydy);    
        
    else
    
        % use nonlinear mass conservation equation (explicit calculation)
        rho0_plus_rho = 2 .* (rhox + rhoy) + rho0;
        rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0_plus_rho .* duxdx);
        rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0_plus_rho .* duydy);
        
    end  
    
    % add in the pre-scaled pressure source term as a mass source    
    if flags.source_p >= t_index
        switch source.p_mode
            case 'dirichlet'
            
                % enforce the source values as a dirichlet boundary condition
                rhox(p_source_pos_index) = source.p(p_source_sig_index, t_index);
                rhoy(p_source_pos_index) = source.p(p_source_sig_index, t_index);
                
            case 'additive'
            
                % extract the source values into a matrix
                source_mat = castZeros(size(rhox));
                source_mat(p_source_pos_index) = source.p(p_source_sig_index, t_index);

                % apply the k-space correction
                switch radial_symmetry
                    case {'WSWA-FFT', 'WSWS-FFT'}
                        
                        % mirror the source matrix in the radial dimension
                        % using the appropriate symmetry 
                        switch radial_symmetry
                            case 'WSWA-FFT'
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         source_mat;
                                source_mat_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2)     = -fliplr(source_mat(:, 2:end));
                                source_mat_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)     =        -source_mat;
                                source_mat_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4)     =  fliplr(source_mat(:, 2:end));
                            case 'WSWS-FFT'
                                source_mat_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         source_mat;
                                source_mat_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2 - 2) =  fliplr(source_mat(:, 2:end - 1));
                        end
                        
                        % apply source correction
                        source_mat_exp = real(ifft2(source_kappa .* fft2(source_mat_exp)));
                        
                        % trim source to correct size
                        source_mat = source_mat_exp(:, 1:kgrid.Ny);
                        
                    case 'WSWA'
                        
                        % compute forward transform in the radial direction
                        % using DCT3 and inverse transform using DCT3^-1 =
                        % DCT2 (WSWA symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DCT3, 2), [], 1)...
                            , [], 1)), DCT2, 2) ./ M;                                         
                        
                    case 'WSWS'
                        
                        % compute forward transform in the radial direction
                        % using DCT1 and inverse transform using DCT1^-1 =
                        % DCT1 (WSWS symmetry), with transforms in the
                        % axial direction computed using the FFT
                        source_mat = dtt1D(real(ifft( ...
                            source_kappa .* fft(dtt1D(source_mat, DCT1, 2), [], 1)...
                            , [], 1)), DCT1, 2) ./ M;                            
                        
                end

                % add the source values to the existing field values, including
                % the k-space correction
                rhox = rhox + source_mat;
                rhoy = rhoy + source_mat;
                
            case 'additive-no-correction'
                
                % add the source values to the existing field values
                rhox(p_source_pos_index) = rhox(p_source_pos_index) + source.p(p_source_sig_index, t_index);
                rhoy(p_source_pos_index) = rhoy(p_source_pos_index) + source.p(p_source_sig_index, t_index);       
            
        end
    end
    
    if ~flags.nonlinear
        switch equation_of_state
            case 'lossless'
                
                % calculate p using a linear adiabatic equation of state
                p = c0.^2 .* (rhox + rhoy);
               
            case 'stokes'
                
                % calculate p using a linear absorbing equation of state
                % assuming alpha_power = 2
                p = c0.^2 .* (...
                    (rhox + rhoy) ...
                    + absorb_tau .* rho0 .* (duxdx + duydy) ...
                    );               
               
        end
    else
        switch equation_of_state
            case 'lossless'
                
                % calculate p using a nonlinear adiabatic equation of state
                p = c0.^2 .* (rhox + rhoy + medium.BonA .* (rhox + rhoy).^2 ./ (2 .* rho0));

            case 'stokes'
                
                % calculate p using a nonlinear absorbing equation of state
                % assuming alpha_power = 2
                p = c0.^2 .* (...
                    (rhox + rhoy) ...
                    + absorb_tau .* rho0 .* (duxdx + duydy) ...
                    + medium.BonA .* (rhox + rhoy).^2 ./ (2 .* rho0) ...
                    );
                
        end
    end
    
    % enforce initial conditions if source.p0 is defined instead of time
    % varying sources
    if t_index == 1 && flags.source_p0
    
        % add the initial pressure to rho as a mass source
        p    = source.p0;
        rhox = source.p0 ./ (2 .* c0.^2);
        rhoy = source.p0 ./ (2 .* c0.^2);
        
        % compute gradient of the pressure
        switch radial_symmetry          
            case {'WSWA-FFT', 'WSWS-FFT'}

                % mirror in radial dimension using appropriate symmetry
                switch radial_symmetry
                    case 'WSWA-FFT'
                        p_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         p;
                        p_exp(:, kgrid.Ny*1 + 2:kgrid.Ny*2)     = -fliplr(p(:, 2:end));
                        p_exp(:, kgrid.Ny*2 + 1:kgrid.Ny*3)     =        -p;
                        p_exp(:, kgrid.Ny*3 + 2:kgrid.Ny*4)     =  fliplr(p(:, 2:end));
                    case 'WSWS-FFT'
                        p_exp(:, kgrid.Ny*0 + 1:kgrid.Ny*1)     =         p;
                        p_exp(:, kgrid.Ny*1 + 1:kgrid.Ny*2 - 2) =  fliplr(p(:, 2:end - 1));
                end

                % compute gradients
                p_k = kappa .* fft2(p_exp);
                dpdx_sgx_exp = real(ifft2( bsxfun(@times, ddx_k_shift_pos, p_k) ));
                dpdy_sgy_exp = real(ifft2( bsxfun(@times, ddy_k .* y_shift_pos, p_k) ));

                % trim gradients
                dpdx_sgx = dpdx_sgx_exp(:, 1:kgrid.Ny);
                dpdy_sgy = dpdy_sgy_exp(:, 1:kgrid.Ny);
                
            case 'WSWA'
                
                % compute forward transform
                p_k = kappa .* fft(dtt1D(p, DCT3, 2), [], 1);

                % compute axial derivative
                dpdx_sgx = dtt1D(real(ifft(bsxfun(@times, ddx_k_shift_pos, p_k), [], 1)), DCT2, 2) ./ M;

                % compute radial derivative 
                dpdy_sgy = dtt1D(real(ifft(bsxfun(@times, ddy_k_wswa, p_k), [], 1)), DST4, 2) ./ M;

            case 'WSWS'

                % compute forward transform
                p_k = kappa .* fft(dtt1D(p, DCT1, 2), [], 1);

                % compute axial derivative
                dpdx_sgx = dtt1D(real(ifft(bsxfun(@times, ddx_k_shift_pos, p_k), [], 1)), DCT1, 2) ./ M;

                % compute radial derivative
                dpdy_sgy(:, 1:end - 1) = dtt1D(real(ifft(bsxfun(@times, ddy_k_wsws, p_k(:, 2:end)), [], 1)), DST3, 2) ./ M;

                % mirror the final point according to HAHA symmetry
                dpdy_sgy(:, end) = -dpdy_sgy(:, end - 1);

        end
        
        % compute u(t = t1 - dt/2) based on u(dt/2) = -u(-dt/2) which
        % forces u(t = t1) = 0         
        ux_sgx = dt .* rho0_sgx_inv .* dpdx_sgx / 2;
        uy_sgy = dt .* rho0_sgy_inv .* dpdy_sgy / 2;

    end  

    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if flags.use_sensor && ~flags.time_rev && (t_index >= sensor.record_start_index)
    
        % update index for data storage 
        file_index = t_index - sensor.record_start_index + 1;
        
        % run sub-function to extract the required data from the acoustic
        % variables
        sensor_data = kspaceFirstOrder_extractSensorData(2, sensor_data, file_index, sensor_mask_index, flags, record, p, ux_sgx, uy_sgy, []);
        
        % extract acoustic pressure with directional response if required -
        % this function needs kgrid and fft2(p), so is not calculated
        % within _extractSensorData
        if flags.compute_directivity && (flags.record_p || flags.record_I || flags.record_I_avg)
            sensor_data.p(:, file_index) = directionalResponse(kgrid, sensor, sensor_mask_index, fft2(p));
        end
        
    end
    
    % estimate the time to run the simulation
    if t_index == ESTIMATE_SIM_TIME_STEPS
        
        % display estimated simulation time
        disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time) * index_end / t_index) '...']);
        
        % check memory usage
        kspaceFirstOrder_checkMemoryUsage;
        
    end    
    
    % plot data if required
    if flags.plot_sim && (rem(t_index, plot_freq) == 0 || t_index == 1 || t_index == index_end) 

        % update progress bar
        waitbar(t_index / kgrid.Nt, pbar);
        drawnow;   

        % ensure p is cast as a CPU variable and remove the PML from the
        % plot if required
        if strcmp(data_cast, 'gpuArray')
            p_plot = double(gather(p(x1:x2, y1:y2)));
        else
            p_plot = double(p(x1:x2, y1:y2));          
        end

        % update plot scale if set to automatic or log
        if flags.plot_scale_auto || flags.plot_scale_log
            kspaceFirstOrder_adjustPlotScale;
        end    

        % add display mask onto plot
        if strcmp(display_mask, 'default')
            p_plot(double(sensor.mask(x1:x2, y1:y2)) == 1) = plot_scale(2);
        elseif ~strcmp(display_mask, 'off')
            p_plot(display_mask(x1:x2, y1:y2) ~= 0) = plot_scale(2);
        end

        if flags.mesh_plot
                       
            % update plot using a mesh without axes
            mesh(y_vec(y1:y2) * scale, kgrid.x_vec(x1:x2) * scale, p_plot, 'EdgeColor', 'Black');
            
            % set aspect ratio and assign plot scale as z limits
            set(gca, 'ZLim', [plot_scale(1)/2, plot_scale(2)]);
            set(gca, 'DataAspectRatio', [1, 1, (plot_scale(2) - plot_scale(1)) / (max(kgrid.x_size, kgrid.y_size) * scale) * 2 ]);
            axis off;
            
            % set view
            view([-38, 40]);
                        
        else
            
            % update plot
            imagesc(y_vec(y1:y2) * scale, kgrid.x_vec(x1:x2) * scale, p_plot, plot_scale);
            colormap(COLOR_MAP);
            ylabel(['axial position [' prefix 'm]']);
            xlabel(['radial position [' prefix 'm]']);
            axis image;
            
        end

        % force plot update
        drawnow;        

        % save movie frame if required
        if flags.record_movie

            % set background color to white
            set(gcf, 'Color', [1, 1, 1]);

            % save the movie frame
            writeVideo(video_obj, getframe(gcf));

        end
        
        % update variable used for timing variable to exclude the first
        % time step if plotting is enabled
        if t_index == 1
            loop_start_time = clock;
        end
        
    end
end

% assign the final time reversal values
if flags.time_rev
    p(sensor_mask_index) = sensor.time_reversal_boundary_data(:, index_end + 1);
end

% update command line status
disp(['  simulation completed in ' scaleTime(toc)]);

% =========================================================================
% CLEAN UP
% =========================================================================

% clean up used figures
if flags.plot_sim
    close(img);
    close(pbar);
    drawnow;
end

% save the movie frames to disk
if flags.record_movie
    close(video_obj);
end

% save the final pressure field if in time reversal mode
if flags.time_rev
    flags.record_p_final = true;
end

% save the final acoustic pressure if required
if flags.record_p_final
    sensor_data.p_final = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
end

% save the final particle velocity if required
if flags.record_u_final
    sensor_data.ux_final = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
    sensor_data.uy_final = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);    
end

% run subscript to cast variables back to double precision if required
if flags.data_recast
    kspaceFirstOrder_dataRecast;
end

% run subscript to compute and save intensity values
if flags.use_sensor && ~flags.time_rev && (flags.record_I || flags.record_I_avg)
    kspaceFirstOrder_saveIntensity;
end

% reorder the sensor points if a binary sensor mask was used for Cartesian
% sensor mask nearest neighbour interpolation (this is performed after
% recasting as the GPU toolboxes do not all support this subscript)
if flags.use_sensor && flags.reorder_data
    kspaceFirstOrder_reorderCartData;
end

% filter the recorded time domain pressure signals if transducer filter
% parameters are given 
if flags.use_sensor && ~flags.time_rev && isfield(sensor, 'frequency_response')
    sensor_data.p = gaussianFilter(sensor_data.p, 1/kgrid.dt, sensor.frequency_response(1), sensor.frequency_response(2));
end

% reorder the sensor points if cuboid corners is used (outputs are indexed
% as [X, Y, T] or [X, Y] rather than [sensor_index, time_index]
if flags.cuboid_corners
    kspaceFirstOrder_reorderCuboidCorners;
end

if ~flags.use_sensor
    
    % if sensor is not used, return empty sensor data
    sensor_data = [];
    
elseif flags.time_rev
    
    % if computing time reversal, reassign sensor_data.p_final to
    % sensor_data
    sensor_data = sensor_data.p_final;
    
elseif ~isfield(sensor, 'record') && ~flags.cuboid_corners
    
    % if sensor.record is not given by the user, and not using a cuboid
    % sensor mask, reassign sensor_data.p to sensor_data
    sensor_data = sensor_data.p;    
    
end

% update command line status
disp(['  total computation time ' scaleTime(etime(clock, start_time))]);

% switch off log
if flags.create_log
    diary off;
end