function [sensor_data, mem_usage] = kspaceFirstOrder3D(kgrid, medium, source, sensor, varargin) %#ok<STOUT>
%KSPACEFIRSTORDER3D 3D time-domain simulation of wave propagation.
%
% DESCRIPTION:
%     kspaceFirstOrder3D simulates the time-domain propagation of
%     compressional waves through a three-dimensional homogeneous or
%     heterogeneous acoustic medium given four input structures: kgrid,
%     medium, source, and sensor. The computation is based on a first-order
%     k-space model which accounts for power law absorption and a
%     heterogeneous sound speed and density. If medium.BonA is specified,
%     cumulative nonlinear effects are also modelled. At each time-step
%     (defined by kgrid.dt and kgrid.Nt or kgrid.t_array), the acoustic
%     field parameters at the positions defined by sensor.mask are recorded
%     and stored. If kgrid.t_array is set to 'auto', this array is
%     automatically generated using the makeTime method of the kWaveGrid
%     class. An anisotropic absorbing boundary layer called a perfectly
%     matched layer (PML) is implemented to prevent waves that leave one
%     side of the domain being reintroduced from the opposite side (a
%     consequence of using the FFT to compute the spatial derivatives in
%     the wave equation). This allows infinite domain simulations to be
%     computed using small computational grids.
%
%     For a homogeneous medium the formulation is exact and the time-steps
%     are only limited by the effectiveness of the perfectly matched layer.
%     For a heterogeneous medium, the solution represents a leap-frog
%     pseudospectral method with a k-space correction that improves the
%     accuracy of computing the temporal derivatives. This allows larger
%     time-steps to be taken for the same level of accuracy compared to
%     conventional pseudospectral time-domain methods. The computational
%     grids are staggered both spatially and temporally.
%
%     An initial pressure distribution can be specified by assigning a
%     matrix (the same size as the computational grid) of arbitrary numeric
%     values to source.p0. A time varying pressure source can similarly be
%     specified by assigning a binary matrix (i.e., a matrix of 1's and 0's
%     with the same dimensions as the computational grid) to source.p_mask
%     where the 1's represent the grid points that form part of the source.
%     The time varying input signals are then assigned to source.p. This
%     can be a single time series (in which case it is applied to all
%     source elements), or a matrix of time series following the source
%     elements using MATLAB's standard column-wise linear matrix index
%     ordering. A time varying velocity source can be specified in an
%     analogous fashion, where the source location is specified by
%     source.u_mask, and the time varying input velocity is assigned to
%     source.ux, source.uy, and source.uz.
%
%     The field values are returned as arrays of time series at the sensor
%     locations defined by sensor.mask. This can be defined in three
%     different ways. (1) As a binary matrix (i.e., a matrix of 1's and 0's
%     with the same dimensions as the computational grid) representing the
%     grid points within the computational grid that will collect the data.
%     (2) As the grid coordinates of two opposing corners of a cuboid in
%     the form [x1; y1; z1; x2; y2; z2]. This is equivalent to using a
%     binary sensor mask covering the same region, however, the output is
%     indexed differently as discussed below. (3) As a series of Cartesian
%     coordinates within the grid which specify the location of the
%     pressure values stored at each time step. If the Cartesian
%     coordinates don't exactly match the coordinates of a grid point, the
%     output values are calculated via interpolation. The Cartesian points
%     must be given as a 3 by N matrix corresponding to the x, y, and z
%     positions, respectively, where the Cartesian origin is assumed to be
%     in the center of the grid. If no output is required, the sensor input
%     can be replaced with an empty array []. Both the source and sensor
%     inputs can also be replaced by an object of the kWaveTransducer
%     class.
%
%     If sensor.mask is given as a set of Cartesian coordinates, the
%     computed sensor_data is returned in the same order. If sensor.mask is
%     given as a binary matrix, sensor_data is returned using MATLAB's
%     standard column-wise linear matrix index ordering. In both cases, the
%     recorded data is indexed as sensor_data(sensor_point_index,
%     time_index). For a binary sensor mask, the field values at a
%     particular time can be restored to the sensor positions within the
%     computation grid using unmaskSensorData. If sensor.mask is given as a
%     list of cuboid corners, the recorded data is indexed as
%     sensor_data(cuboid_index).p(x_index, y_index, z_index, time_index),
%     where x_index, y_index, and z_index correspond to the grid index
%     within the cuboid, and cuboid_index corresponds to the number of the
%     cuboid if more than one is specified.
%
%     By default, the recorded acoustic pressure field is passed directly
%     to the output sensor_data. However, other acoustic parameters can
%     also be recorded by setting sensor.record to a cell array of the form
%     {'p', 'u', 'p_max', ...}. For example, both the particle velocity and
%     the acoustic pressure can be returned by setting sensor.record =
%     {'p', 'u'}. If sensor.record is given, the output sensor_data is
%     returned as a structure with the different outputs appended as
%     structure fields. For example, if sensor.record = {'p', 'p_final',
%     'p_max', 'u'}, the output would contain fields sensor_data.p,
%     sensor_data.p_final, sensor_data.p_max, sensor_data.ux,
%     sensor_data.uy, and sensor_data.uz. Most of the output parameters are
%     recorded at the given sensor positions and are indexed as
%     sensor_data.field(sensor_point_index, time_index) or
%     sensor_data(cuboid_index).field(x_index, y_index, z_index,
%     time_index) if using a sensor mask defined as cuboid corners. The
%     exceptions are the averaged quantities ('p_max', 'p_rms', 'u_max',
%     'p_rms', 'I_avg'), the 'all' quantities ('p_max_all', 'p_min_all',
%     'u_max_all', 'u_min_all'), and the final quantities ('p_final',
%     'u_final'). The averaged quantities are indexed as
%     sensor_data.p_max(sensor_point_index) or
%     sensor_data(cuboid_index).p_max(x_index, y_index, z_index) if using
%     cuboid corners, while the final and 'all' quantities are returned
%     over the entire grid and are always indexed as
%     sensor_data.p_final(nx, ny, nz), regardless of the type of sensor
%     mask.
%
%     kspaceFirstOrder3D may also be used for time reversal image
%     reconstruction by assigning the time varying pressure recorded over
%     an arbitrary sensor surface to the input field
%     sensor.time_reversal_boundary_data. This data is then enforced in
%     time reversed order as a time varying Dirichlet boundary condition
%     over the sensor surface given by sensor.mask. The boundary data must
%     be indexed as sensor.time_reversal_boundary_data(sensor_point_index,
%     time_index). If sensor.mask is given as a set of Cartesian
%     coordinates, the boundary data must be given in the same order. An
%     equivalent binary sensor mask (computed using nearest neighbour
%     interpolation) is then used to place the pressure values into the
%     computational grid at each time step. If sensor.mask is given as a
%     binary matrix of sensor points, the boundary data must be ordered
%     using MATLAB's standard column-wise linear matrix indexing. If no
%     additional inputs are required, the source input can be replaced with
%     an empty array [].
%
%     Acoustic attenuation compensation can also be included during time
%     reversal image reconstruction by assigning the absorption parameters
%     medium.alpha_coeff and medium.alpha_power and reversing the sign of
%     the absorption term by setting medium.alpha_sign = [-1, 1]. This
%     forces the propagating waves to grow according to the absorption
%     parameters instead of decay. The reconstruction should then be
%     regularised by assigning a filter to medium.alpha_filter (this can be
%     created using getAlphaFilter).
%
%     Note: To run a simple photoacoustic image reconstruction example
%     using time reversal (that commits the 'inverse crime' of using the
%     same numerical parameters and model for data simulation and image
%     reconstruction), the sensor_data returned from a k-Wave simulation
%     can be passed directly to sensor.time_reversal_boundary_data with the
%     input fields source.p0 and source.p removed or set to zero.
%         
% USAGE:
%     sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor)
%     sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
% The minimum fields that must be assigned to run an initial value problem
% (for example, a photoacoustic forward simulation) are marked with a *. 
%
%     kgrid*                 - k-Wave grid object returned by kWaveGrid
%                              containing Cartesian and k-space grid fields 
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
%     medium.alpha_power     - power law absorption exponent
%     medium.alpha_coeff     - power law absorption coefficient 
%                              [dB/(MHz^y cm)] 
%     medium.alpha_mode      - optional input to force either the
%                              absorption or dispersion terms in the
%                              equation of state to be excluded; valid
%                              inputs are 'no_absorption' or
%                              'no_dispersion'  
%     medium.alpha_filter    - frequency domain filter applied to the
%                              absorption and dispersion terms in the
%                              equation of state 
%     medium.alpha_sign      - two element array used to control the sign
%                              of absorption and dispersion terms in the
%                              equation of state  
%
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
%                              x-direction at each of the source positions
%                              given by source.u_mask
%     source.uy              - time varying particle velocity in the
%                              y-direction at each of the source positions
%                              given by source.u_mask
%     source.uz              - time varying particle velocity in the
%                              z-direction at each of the source positions
%                              given by source.u_mask  
%     source.u_mask          - binary matrix specifying the positions of
%                              the time varying particle velocity
%                              distribution
%     source.u_mode          - optional input to control whether the input
%                              velocity is applied as a force source or
%                              enforced as a dirichlet boundary condition;
%                              valid inputs are 'additive' (the default) or
%                              'dirichlet'
%
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
%                              kspaceFirstOrder3DC. Can be set to an
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
%     'MovieArgs'            - Settings for VideoWriter. Parameters must be
%                              given as {'param', value, ...} pairs within
%                              a cell array (default = {}), where 'param'
%                              corresponds to a writable property of a
%                              VideoWriter object. 
%     'MovieName'            - Name of the movie produced when
%                              'RecordMovie' is set to true (default =
%                              'date-time-kspaceFirstOrder3D'). 
%     'MovieProfile'         - Profile input passed to VideoWriter.
%     'PlotFreq'             - The number of iterations which must pass 
%                              before the simulation plot is updated
%                              (default = 10).
%     'PlotLayout'           - Boolean controlling whether three plots are
%                              produced of the initial simulation layout
%                              (initial pressure, sound speed, density)
%                              (default = false). 
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
%                              PMLSize and PMLAlpha can be given as three
%                              element arrays to specify the x, y, and z
%                              properties, respectively. To remove the PML,
%                              set the appropriate PMLAlpha to zero rather
%                              than forcing the PML to be of zero size
%                              (default = 10).
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
%     'SaveToDisk'           - String containing a filename (including 
%                              pathname if required). If set, after the
%                              precomputation phase, the input variables
%                              used in the time loop are saved the
%                              specified location in HDF5 format. The
%                              simulation then exits. The saved variables
%                              can be used to run simulations using the C++
%                              code.
%     'StreamToDisk'         - Boolean controlling whether sensor_data is
%                              periodically saved to disk to avoid storing
%                              the complete matrix in memory. StreamToDisk
%                              may also be given as an integer which
%                              specifies the number of times steps that are
%                              taken before the data is saved to disk
%                              (default = 200).
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
%     sensor_data.uz         - time varying particle velocity in the
%                              z-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'u' is
%                              set)
%     sensor_data.ux_max     - maximum particle velocity in the x-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_max' is set)     
%     sensor_data.uy_max     - maximum particle velocity in the y-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_max' is set) 
%     sensor_data.uz_max     - maximum particle velocity in the z-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_max' is set) 
%     sensor_data.ux_min     - minimum particle velocity in the x-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_min' is set)
%     sensor_data.uy_min     - minimum particle velocity in the y-direction
%                              recorded at the sensor positions given by
%                              sensor.mask (returned if 'u_min' is set)
%     sensor_data.uz_min     - minimum particle velocity in the z-direction
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
%     sensor_data.uz_rms     - rms of the time varying particle velocity in
%                              the z-direction recorded at the sensor
%                              positions given by sensor.mask (returned if
%                              'u_rms' is set) 
%     sensor_data.ux_final   - final particle velocity field in the
%                              x-direction at all grid points within the
%                              domain (returned if 'u_final' is set) 
%     sensor_data.uy_final   - final particle velocity field in the
%                              y-direction at all grid points within the
%                              domain (returned if 'u_final' is set)  
%     sensor_data.uz_final   - final particle velocity field in the
%                              z-direction at all grid points within the
%                              domain (returned if 'u_final' is set) 
%     sensor_data.ux_max_all - maximum particle velocity in the x-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_max_all' is set) 
%     sensor_data.uy_max_all - maximum particle velocity in the y-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_max_all' is set) 
%     sensor_data.uz_max_all - maximum particle velocity in the z-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_max_all' is set) 
%     sensor_data.ux_min_all - minimum particle velocity in the x-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_min_all' is set)
%     sensor_data.uy_min_all - minimum particle velocity in the y-direction
%                              recorded at all grid points within the
%                              domain (returned if 'u_min_all' is set) 
%     sensor_data.uz_min_all - minimum particle velocity in the z-direction
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
%     sensor_data.uz_non_staggered 
%                            - time varying particle velocity in the
%                              z-direction recorded at the sensor positions
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
%     sensor_data.Iz         - time varying acoustic intensity in the
%                              z-direction recorded at the sensor positions
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
%     sensor_data.Iz_avg     - average acoustic intensity in the
%                              z-direction recorded at the sensor positions
%                              given by sensor.mask (returned if 'I_avg' is
%                              set)  
%
% ABOUT:
%     author                 - Bradley Treeby and Ben Cox
%     date                   - 7th April 2009
%     last update            - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby and Ben Cox
%
% See also kspaceFirstOrderAS, kspaceFirstOrder1D, kspaceFirstOrder2D,
% kspaceFirstOrder3DC, kspaceFirstOrder3DG, kWaveGrid, kWaveTransducer,
% kspaceSecondOrder

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
% where sgx = (x + dx/2, y, z), sgy = (x, y + dy/2, z), sgz = (x, y, z +
% dz/2)
if numDim(rho0) == 3 && flags.use_sg
    
    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x + kgrid.dx/2, kgrid.y, kgrid.z, '*linear');
    rho0_sgy = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y + kgrid.dy/2, kgrid.z, '*linear');
    rho0_sgz = interpn(kgrid.x, kgrid.y, kgrid.z, rho0, kgrid.x, kgrid.y, kgrid.z + kgrid.dz/2, '*linear');
    
    % set values outside of the interpolation range to original values
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx));
    rho0_sgy(isnan(rho0_sgy)) = rho0(isnan(rho0_sgy));    
    rho0_sgz(isnan(rho0_sgz)) = rho0(isnan(rho0_sgz));
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
    rho0_sgy = rho0;
    rho0_sgz = rho0;
end

% invert rho0 so it doesn't have to be done each time step
rho0_sgx_inv = 1./rho0_sgx;
rho0_sgy_inv = 1./rho0_sgy;
rho0_sgz_inv = 1./rho0_sgz;

% clear unused variables if not using them in _saveToDisk
if ~flags.save_to_disk
    clear rho0_sgx rho0_sgy rho0_sgz
end

% =========================================================================
% PREPARE DERIVATIVE AND PML OPERATORS
% =========================================================================

% get the PML operators based on the reference sound speed and PML settings
pml_x     = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, false, 1);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, true && flags.use_sg, 1);
pml_y     = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, false, 2);
pml_y_sgy = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, true && flags.use_sg, 2);
pml_z     = getPML(kgrid.Nz, kgrid.dz, kgrid.dt, c_ref, pml_z_size, pml_z_alpha, false, 3);
pml_z_sgz = getPML(kgrid.Nz, kgrid.dz, kgrid.dt, c_ref, pml_z_size, pml_z_alpha, true && flags.use_sg, 3);

% define the k-space derivative operators, multiply by the staggered
% grid shift operators, and then re-order using ifftshift (the option
% flags.use_sg exists for debugging) 
if flags.use_sg
    ddx_k_shift_pos = ifftshift( 1i * kgrid.kx_vec .* exp( 1i * kgrid.kx_vec * kgrid.dx/2) );
    ddx_k_shift_neg = ifftshift( 1i * kgrid.kx_vec .* exp(-1i * kgrid.kx_vec * kgrid.dx/2) );
    ddy_k_shift_pos = ifftshift( 1i * kgrid.ky_vec .* exp( 1i * kgrid.ky_vec * kgrid.dy/2) );
    ddy_k_shift_neg = ifftshift( 1i * kgrid.ky_vec .* exp(-1i * kgrid.ky_vec * kgrid.dy/2) );
    ddz_k_shift_pos = ifftshift( 1i * kgrid.kz_vec .* exp( 1i * kgrid.kz_vec * kgrid.dz/2) );
    ddz_k_shift_neg = ifftshift( 1i * kgrid.kz_vec .* exp(-1i * kgrid.kz_vec * kgrid.dz/2) );
else
    ddx_k_shift_pos = ifftshift( 1i * kgrid.kx_vec );
    ddx_k_shift_neg = ifftshift( 1i * kgrid.kx_vec );
    ddy_k_shift_pos = ifftshift( 1i * kgrid.ky_vec );
    ddy_k_shift_neg = ifftshift( 1i * kgrid.ky_vec );
    ddz_k_shift_pos = ifftshift( 1i * kgrid.kz_vec );
    ddz_k_shift_neg = ifftshift( 1i * kgrid.kz_vec );         
end

% force the derivative and shift operators to be in the correct direction
% for use with BSXFUN
ddy_k_shift_pos = ddy_k_shift_pos.'; 
ddy_k_shift_neg = ddy_k_shift_neg.';
ddz_k_shift_pos = permute(ddz_k_shift_pos, [2, 3, 1]);
ddz_k_shift_neg = permute(ddz_k_shift_neg, [2, 3, 1]);

% create k-space operators (the option flags.use_kspace exists for debugging)
if flags.use_kspace
    kappa = ifftshift(sinc(c_ref .* kgrid.k .* kgrid.dt / 2));
    if (flags.source_p && strcmp(source.p_mode, 'additive')) || ((flags.source_ux || flags.source_uy || flags.source_uz) && strcmp(source.u_mode, 'additive'))
        source_kappa = ifftshift(cos(c_ref .* kgrid.k .* kgrid.dt / 2));
    end
else
    kappa        = 1;
    source_kappa = 1;
end

% =========================================================================
% SAVE DATA TO DISK FOR RUNNING SIMULATION EXTERNAL TO MATLAB
% =========================================================================

% save to disk option for saving the input matrices to disk for running
% simulations using k-Wave++
if flags.save_to_disk
    
    % store the pml size for resizing transducer object below
    retract_size = [pml_x_size, pml_y_size, pml_z_size]; 
    
    % run subscript to save files to disk
    kspaceFirstOrder_saveToDisk;
    
    % run subscript to resize the transducer object if the grid has been
    % expanded 
    kspaceFirstOrder_retractTransducerGridSize;
    
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
p      = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhox   = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhoy   = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
rhoz   = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
ux_sgx = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
uy_sgy = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
uz_sgz = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
p_k    = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);

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
    index_step = 1;
    index_end = kgrid.Nt; 
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
    [x_sc, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec; kgrid.z_vec]));  %#ok<ASGLU>
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

        % update p_k
        p_k = fftn(p);

        % compute rhox and rhoz using an adiabatic equation of state
        rhox_mod = p ./ (3 .* c0.^2);
        rhoy_mod = p ./ (3 .* c0.^2);
        rhoz_mod = p ./ (3 .* c0.^2);
        rhox(sensor_mask_index) = rhox_mod(sensor_mask_index);
        rhoy(sensor_mask_index) = rhoy_mod(sensor_mask_index);
        rhoz(sensor_mask_index) = rhoz_mod(sensor_mask_index);
           
    end    

    % calculate ux, uy and uz at the next time step using dp/dx, dp/dy and
    % dp/dz at the current time step
    ux_sgx = bsxfun(@times, pml_x_sgx, ...
        bsxfun(@times, pml_x_sgx, ux_sgx) ... 
        - dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* p_k) )) ...
        );
    uy_sgy = bsxfun(@times, pml_y_sgy, ...
        bsxfun(@times, pml_y_sgy, uy_sgy) ...
        - dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* p_k) )) ...
        );
    uz_sgz = bsxfun(@times, pml_z_sgz, ...
        bsxfun(@times, pml_z_sgz, uz_sgz) ...
        - dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* p_k) )) ...
        );                   
    
    % add in the velocity source terms   
    if flags.source_ux >= t_index
        switch source.u_mode
            case 'dirichlet'
            
                % enforce the source values as a dirichlet boundary condition
                ux_sgx(u_source_pos_index) = source.ux(u_source_sig_index, t_index);   
                
            case 'additive'
                
                % extract the source values into a matrix
                source_mat = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
                source_mat(u_source_pos_index) = source.ux(u_source_sig_index, t_index);

                % apply the k-space correction
                source_mat = real(ifftn(source_kappa .* fftn(source_mat)));

                % add the source values to the existing field values
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
                source_mat = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
                source_mat(u_source_pos_index) = source.uy(u_source_sig_index, t_index);

                % apply the k-space correction
                source_mat = real(ifftn(source_kappa .* fftn(source_mat)));

                % add the source values to the existing field values
                % including the k-space correction
                uy_sgy = uy_sgy + source_mat;                    
                
            case 'additive-no-correction'

                % add the source values to the existing field values
                uy_sgy(u_source_pos_index) = uy_sgy(u_source_pos_index) + source.uy(u_source_sig_index, t_index);
            
        end        
    end
    if flags.source_uz >= t_index
        switch source.u_mode
            case 'dirichlet'  
            
                % enforce the source values as a dirichlet boundary condition
                uz_sgz(u_source_pos_index) = source.uz(u_source_sig_index, t_index);  
            
            case 'additive'                
                
                % extract the source values into a matrix
                source_mat = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
                source_mat(u_source_pos_index) = source.uz(u_source_sig_index, t_index);

                % apply the k-space correction
                source_mat = real(ifftn(source_kappa .* fftn(source_mat)));

                % add the source values to the existing field values
                % including the k-space correction
                uz_sgz = uz_sgz + source_mat;                    
                
            case 'additive-no-correction'
            
                % add the source values to the existing field values
                uz_sgz(u_source_pos_index) = uz_sgz(u_source_pos_index) + source.uz(u_source_sig_index, t_index);
            
        end         
    end   
    
    % add in transducer source term; there will normally be less source
    % points than time points, so this is only done until the source points
    % run out (the number of source points is stored in flags.transducer_source)
    if flags.transducer_source >= t_index
        
        % as only flat transducers are currently supported, assume all the
        % energy is transfered to x-direction velocity, multiply source
        % terms by apodization weights
        ux_sgx(u_source_pos_index) = ux_sgx(u_source_pos_index) + transducer_transmit_apodization.*transducer_input_signal(delay_mask);
        
        % update the delay_mask - this maps the source positions belonging
        % to different transducer elements to time points within
        % transducer_input_signal (this is a single time series) based on the
        % beamforming and focussing delays.  
        delay_mask = delay_mask + 1;
        
    end
    
    % calculate dux/dx, duydy and duz/dz at the next time step
    duxdx = real(ifftn( bsxfun(@times, ddx_k_shift_neg, kappa .* fftn(ux_sgx)) ));
    duydy = real(ifftn( bsxfun(@times, ddy_k_shift_neg, kappa .* fftn(uy_sgy)) ));
    duzdz = real(ifftn( bsxfun(@times, ddz_k_shift_neg, kappa .* fftn(uz_sgz)) ));        

    % calculate rhox, rhoy and rhoz at the next time step
    if ~flags.nonlinear
        
        % use linearised mass conservation equation
        rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0 .* duxdx);
        rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0 .* duydy);        
        rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0 .* duzdz);
        
    else

        % use nonlinear mass conservation equation (explicit calculation)
        rho0_plus_rho = 2 .* (rhox + rhoy + rhoz) + rho0;
        rhox = bsxfun(@times, pml_x, bsxfun(@times, pml_x, rhox) - dt .* rho0_plus_rho .* duxdx);
        rhoy = bsxfun(@times, pml_y, bsxfun(@times, pml_y, rhoy) - dt .* rho0_plus_rho .* duydy);
        rhoz = bsxfun(@times, pml_z, bsxfun(@times, pml_z, rhoz) - dt .* rho0_plus_rho .* duzdz);
                
    end 
    
    % add in the pre-scaled pressure source term as a mass source   
    if flags.source_p >= t_index
        switch source.p_mode
            case 'dirichlet'
            
                % enforce source values as a dirichlet boundary condition
                rhox(p_source_pos_index) = source.p(p_source_sig_index, t_index);  
                rhoy(p_source_pos_index) = source.p(p_source_sig_index, t_index);  
                rhoz(p_source_pos_index) = source.p(p_source_sig_index, t_index);
            
            case 'additive'

                % extract the source values into a matrix
                source_mat = castZeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
                source_mat(p_source_pos_index) = source.p(p_source_sig_index, t_index);

                % apply the k-space correction
                source_mat = real(ifftn(source_kappa .* fftn(source_mat)));

                % add the source values to the existing field values
                % including the k-space correction
                rhox = rhox + source_mat;
                rhoy = rhoy + source_mat;
                rhoz = rhoz + source_mat;
                
            case 'additive-no-correction'
            
                % add the source values to the existing field values
                rhox(p_source_pos_index) = rhox(p_source_pos_index) + source.p(p_source_sig_index, t_index);
                rhoy(p_source_pos_index) = rhoy(p_source_pos_index) + source.p(p_source_sig_index, t_index);  
                rhoz(p_source_pos_index) = rhoz(p_source_pos_index) + source.p(p_source_sig_index, t_index);
            
        end
    end
    
    % calculate p at the next time step
    if ~flags.nonlinear
        switch equation_of_state
            case 'lossless'
                
                % calculate p using a linear adiabatic equation of state
                p = c0.^2 .* (rhox + rhoy + rhoz);
                
            case 'absorbing'
                
                % calculate p using a linear absorbing equation of state                
                p = c0.^2 .* ( ...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) )) ...
                    - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) )) ...
                    );  
               
            case 'stokes'
                
                % calculate p using a nonlinear absorbing equation of state
                % assuming alpha_power = 2
                p = c0.^2 .* (...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau .* rho0 .* (duxdx + duydy + duzdz) ...
                    );
               
        end
    else
        switch equation_of_state
            case 'lossless'
                
                % calculate p using a nonlinear adiabatic equation of state
                p = c0.^2 .* (rhox + rhoy + rhoz + medium.BonA .* (rhox + rhoy + rhoz).^2 ./ (2 .* rho0));
                
            case 'absorbing'
                
                % calculate p using a nonlinear absorbing equation of state 
                p = c0.^2 .* (...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau .* real(ifftn( absorb_nabla1 .* fftn(rho0 .* (duxdx + duydy + duzdz)) ))...
                    - absorb_eta .* real(ifftn( absorb_nabla2 .* fftn(rhox + rhoy + rhoz) ))...
                    + medium.BonA .*(rhox + rhoy + rhoz).^2 ./ (2 .* rho0) ...
                    );

            case 'stokes'
                
                % calculate p using a nonlinear absorbing equation of state
                % assuming alpha_power = 2
                p = c0.^2 .* (...
                    (rhox + rhoy + rhoz) ...
                    + absorb_tau .* rho0 .* (duxdx + duydy + duzdz) ...
                    + medium.BonA .* (rhox + rhoy + rhoz).^2 ./ (2 .* rho0) ...
                    ); 
                
        end
    end
             
    % enforce initial conditions if source.p0 is defined instead of time
    % varying sources
    if t_index == 1 && flags.source_p0
    
        % add the initial pressure to rho as a mass source
        p = source.p0;
        rhox = source.p0 ./ (3 .* c0.^2);
        rhoy = source.p0 ./ (3 .* c0.^2);
        rhoz = source.p0 ./ (3 .* c0.^2);

        % compute u(t = t1 + dt/2) based on the assumption u(dt/2) = -u(-dt/2)
        % which forces u(t = t1) = 0
        ux_sgx = dt .* rho0_sgx_inv .* real(ifftn( bsxfun(@times, ddx_k_shift_pos, kappa .* fftn(p)) )) / 2;
        uy_sgy = dt .* rho0_sgy_inv .* real(ifftn( bsxfun(@times, ddy_k_shift_pos, kappa .* fftn(p)) )) / 2;
        uz_sgz = dt .* rho0_sgz_inv .* real(ifftn( bsxfun(@times, ddz_k_shift_pos, kappa .* fftn(p)) )) / 2;    
        
    end
    
    % precompute fft of p here so p can be modified for visualisation
    p_k = fftn(p);        

    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if flags.use_sensor && ~flags.time_rev && (t_index >= sensor.record_start_index)
    
        % update index for data storage - if streaming to disk, a smaller
        % matrix is used which is continually overwritten, and then saved to
        % disk each time it is filled
        if flags.stream_to_disk
            file_index = t_index - flags.stream_to_disk * (stream_data_index - 1) - sensor.record_start_index + 1;
        else
            file_index = t_index - sensor.record_start_index + 1;
        end
        
        % store the acoustic pressure if using a transducer object
        if flags.transducer_sensor
            
            % check if an elevation focus is set
            if flags.transducer_receive_elevation_focus

                % update the sensor data buffer
                sensor_data_buffer = circshift(sensor_data_buffer, [0, 1]);
                sensor_data_buffer(:, 1) = p(sensor_mask_index);

                % if buffer has been filled, store average pressure
                % across each transducer element accounting for the
                % elevation focus
                if file_index >= sensor_data_buffer_size
                    
                    % get the current values
                    current_vals = sum(transducer_receive_mask .* sensor_data_buffer, 2);

                    % reshape and average
                    sensor_data.transducer(:, file_index - sensor_data_buffer_size + 1) = ...
                        sum(reshape(sum(reshape(current_vals, sensor.number_active_elements * sensor.element_width, sensor.element_length), 2), sensor.element_width, sensor.number_active_elements).', 2);  
                    
                end
                
            else
                
                % store average pressure across each transducer element
                sensor_data.transducer(:, file_index) = ...
                    sum(reshape(sum(reshape(p(sensor_mask_index), sensor.number_active_elements * sensor.element_width, sensor.element_length), 2), sensor.element_width, sensor.number_active_elements).', 2);
                
            end
        
        % extract data from specified grid points
        else
           
            % run sub-function to extract the required data from the
            % acoustic variables
            sensor_data = kspaceFirstOrder_extractSensorData(3, sensor_data, file_index, sensor_mask_index, flags, record, p, ux_sgx, uy_sgy, uz_sgz);
            
        end
            
        % if the data is being streamed to disk and sensor_data has just been
        % filled, append the values of sensor_data to the values already saved
        % to disk (this option is currently only supported for recording
        % the time varying pressure)
        if flags.stream_to_disk && (file_index == flags.stream_to_disk)

            % open the file to append the new data
            if stream_data_index == 1
                
                % create or open the file and overwrite any existing data
                try
                    fid = fopen(STREAM_TO_DISK_FILENAME, 'w+');
                catch ME
                    disp('Error in writing file using ''StreamToDisk''');
                    rethrow(ME);
                end
                
            else
                
                % open the file and append new data
                fid = fopen(STREAM_TO_DISK_FILENAME, 'a+');
                
            end

            % write values at end of file using the precision specified by
            % the data_cast option. If using the GPU, the PCT requires a
            % call to gather, while GPUmat and Accelereyes only require a
            % cast back to single or double.
            if strcmp(data_cast, 'gpuArray')
                if strcmp(precision, 'single')
                    fwrite(fid, gather(sensor_data.p), 'single');
                else
                    fwrite(fid, gather(sensor_data.p), 'double');
                end
            elseif strcmp(precision, 'single')
                fwrite(fid, single(sensor_data.p), 'single');
            else
                fwrite(fid, double(sensor_data.p), 'double');       
            end

            % close the file
            fclose(fid);

            % increment the data file index if there is
            % still data left
            if index_end > t_index
                stream_data_index = stream_data_index + 1;
            end

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
            p_plot = double(gather(p(x1:x2, y1:y2, z1:z2)));
        else
            p_plot = double(p(x1:x2, y1:y2, z1:z2));       
        end

        % update plot scale if set to automatic or log
        if flags.plot_scale_auto || flags.plot_scale_log
            kspaceFirstOrder_adjustPlotScale;
        end         
        
        % add display mask onto plot
        if strcmp(display_mask, 'default')
            p_plot(sensor.mask(x1:x2, y1:y2, z1:z2) ~= 0) = plot_scale(2);
        elseif ~strcmp(display_mask, 'off')
            p_plot(display_mask(x1:x2, y1:y2, z1:z2) ~= 0) = plot_scale(2);
        end     

        % update plot
        planeplot(scale * kgrid.x_vec(x1:x2), scale * kgrid.y_vec(y1:y2), scale * kgrid.z_vec(z1:z2), p_plot, '', plot_scale, prefix, COLOR_MAP);

        % save movie frames if required
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

% save final sensor_data variable to disk if required
if flags.stream_to_disk

    % clear some time loop variables before reloading sensor data to
    % free up memory (in case the sensor data is very large)
    clear duxdx duydy duzdz ux_sgx uy_sgy uz_sgz rhox rhoy rhoz;
    
    % double check the data has not just been saved
    if file_index ~= flags.stream_to_disk
        
        % open the file to append the new data
        if stream_data_index == 1
            
            % open the file and overwrite any existing data
            fid = fopen(STREAM_TO_DISK_FILENAME, 'w+');
            
        else
            
            % open the file and append new data
            fid = fopen(STREAM_TO_DISK_FILENAME, 'a+');
            
        end
        
        % extract required data
        sensor_data.p = sensor_data.p(:, 1:file_index);
        
        % write values at end of file
        if strcmp(data_cast, 'gpuArray')
            if strcmp(precision, 'single')
                fwrite(fid, gather(sensor_data.p), 'single');
            else
                fwrite(fid, gather(sensor_data.p), 'double');
            end
        elseif strcmp(precision, 'single')
            fwrite(fid, single(sensor_data.p), 'single');
        else
            fwrite(fid, double(sensor_data.p), 'double');       
        end
            
        % close the file
        fclose(fid);
        
    end

    % compute the number of recorded time points
    num_stream_time_points = kgrid.Nt - sensor.record_start_index + 1;
    
    % reload complete streamed data and assign to sensor_data
    fid = fopen(STREAM_TO_DISK_FILENAME, 'r');
    if strcmp(precision, 'single')
        sensor_data.p = fread(fid, [num_sensor_points, num_stream_time_points], '*single');            
    else
        sensor_data.p = fread(fid, [num_sensor_points, num_stream_time_points], '*double');
    end
    fclose(fid);    
    
    % permanently delete the temporary storage
    disp('  removing temporary data...');
    old_state = recycle;
    recycle('off');
    delete(STREAM_TO_DISK_FILENAME);
    recycle(old_state);
   
end

% save the final pressure field if in time reversal mode
if flags.time_rev
    flags.record_p_final = true;
end

% save the final acoustic pressure if required
if flags.record_p_final
    sensor_data.p_final = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
end

% save the final particle velocity if required
if flags.record_u_final
    sensor_data.ux_final = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
    sensor_data.uy_final = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
    sensor_data.uz_final = uz_sgz(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
end

% process the sensor data if recorded using a transducer
if flags.transducer_sensor
    
    % the pressure values across the grid points in each element have
    % already been automatically summed, so now scale the summed pressure
    % data by the number of grid points in each element to convert to an
    % average 
    sensor_data.transducer = sensor_data.transducer ./ (sensor.element_length * sensor.element_width);
    
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
% as [X, Y, Z, T] or [X, Y, Z] rather than [sensor_index, time_index]
if flags.cuboid_corners
    kspaceFirstOrder_reorderCuboidCorners;
end

% run subscript to resize the transducer object if the grid has been
% expanded 
kspaceFirstOrder_retractTransducerGridSize;

if ~flags.use_sensor
    
    % if sensor is not used, return empty sensor data
    sensor_data = [];
    
elseif flags.transducer_sensor
    
    % if using a transducer as sensor, reassign sensor_data.transducer to
    % sensor_data
    sensor_data = sensor_data.transducer;
    
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

function planeplot(x_vec, y_vec, z_vec, data, data_title, plot_scale, prefix, color_map)
% Subfunction to produce a plot of a three-dimensional matrix through the
% three central planes.
   
subplot(2, 2, 1);
imagesc(y_vec, x_vec, squeeze(data(:, :, round(end/2))), plot_scale);
title([data_title 'x-y plane']);
axis image;

subplot(2, 2, 2);
imagesc(z_vec, x_vec, squeeze(data(:, round(end/2), :)), plot_scale);
title('x-z plane');
axis image;
xlabel(['(All axes in ' prefix 'm)']);

subplot(2, 2, 3);
imagesc(z_vec, y_vec, squeeze(data(round(end/2), :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(color_map); 
drawnow;