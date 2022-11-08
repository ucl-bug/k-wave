% DESCRIPTION:
%     Subscript to set the default values of the flags, literals, and
%     settings used within the fluid and elastic simulation codes.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 15th May 2018
%     last update - 18th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby

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

% =========================================================================
% FLAGS WHICH DEPEND ON USER INPUTS (THESE SHOULD NOT BE MODIFIED)
% =========================================================================

% flags which control the type of simulation
flags.elastic_code              = false;                % true if elastic simulation
flags.kspace_elastic_code       = false;                % true if elastic simulation with k-space correction
flags.kelvin_voigt_model        = false;                % true if elastic simulation with absorption
flags.nonuniform_grid           = false;                % true if the computational grid is non-uniform
flags.axisymmetric              = false;                % true if fluid axisymmetric simulation
flags.nonlinear                 = false;                % true if fluid simulation with nonlinearity
flags.absorbing                 = false;                % true if simulation with absorption
flags.stokes                    = false;                % true if simulation with absorption and y = 2

% flags which control the characteristics of the sensor
flags.use_sensor                = false;                % false if no output of any kind is required
flags.blank_sensor              = false;                % true if sensor.mask is not defined but _max_all or _final variables are still recorded
flags.time_rev                  = false;                % true for time reversal simulaions using sensor.time_reversal_boundary_data
flags.elastic_time_rev          = false;                % true if using time reversal with the elastic code
flags.compute_directivity       = false;                % true if directivity calculations in 2D are used by setting sensor.directivity_angle
flags.reorder_data              = false;                % true if sensor.mask is Cartesian with nearest neighbour interpolation which is calculated using a binary mask and thus must be re-ordered
flags.binary_sensor_mask        = true;                 % true if sensor.mask is a binary mask
flags.cuboid_corners            = false;                % true if sensor.mask is a list of cuboid corners
flags.transducer_sensor         = false;                % true if sensor is an object of the kWaveTransducer class

% flags which control which parameters are recorded
flags.record_p                  = true;                 % time-varying pressure
flags.record_p_max              = false;                % maximum pressure over simulation
flags.record_p_min              = false;                % minimum pressure over simulation
flags.record_p_rms              = false;                % root-mean-squared pressure over simulation
flags.record_p_max_all          = false;                % maximum pressure over simulation at all grid points
flags.record_p_min_all          = false;                % minimum pressure over simulation at all grid points
flags.record_p_final            = false;                % final pressure field at all grid points
flags.record_u                  = false;                % time-varying particle velocity
flags.record_u_split_field      = false;                % compressional and shear components of time-varying particle velocity
flags.record_u_non_staggered    = false;                % time-varying particle velocity on non-staggered grid
flags.record_u_max              = false;                % maximum particle velocity over simulation
flags.record_u_min              = false;                % minimum particle velocity over simulation
flags.record_u_rms              = false;                % root-mean-squared particle velocity over simulation
flags.record_u_max_all          = false;                % maximum particle velocity over simulation at all grid points
flags.record_u_min_all          = false;                % minimum particle velocity over simulation at all grid points
flags.record_u_final            = false;                % final particle velocity field at all grid points
flags.record_I                  = false;                % time-varying acoustic intensity
flags.record_I_avg              = false;                % time-averaged acoustic intensity

% flags which control the types of source used
flags.source_p0                 = false;                % initial pressure
flags.source_p0_elastic         = false;                % initial pressure in the elastic code
flags.source_p                  = false;                % time-varying pressure
flags.source_p_labelled         = false;                % time-varying pressure with labelled source mask
flags.source_ux                 = false;                % time-varying particle velocity
flags.source_uy                 = false;                % time-varying particle velocity
flags.source_uz                 = false;                % time-varying particle velocity
flags.source_u_labelled         = false;                % time-varying particle velocity with labelled source mask
flags.source_sxx                = false;                % time-varying stress
flags.source_syy                = false;                % time-varying stress
flags.source_szz                = false;                % time-varying stress
flags.source_sxy                = false;                % time-varying stress
flags.source_sxz                = false;                % time-varying stress
flags.source_syz                = false;                % time-varying stress
flags.source_s_labelled         = false;                % time-varying stress with labelled source mask
flags.use_w_source_correction_p = false;                % use the w source correction instead of the k-space source correction for p sources
flags.use_w_source_correction_u = false;                % use the w source correction instead of the k-space source correction for u sources

% transducer source flags
flags.transducer_source         = false;                % transducer is object of kWaveTransducer class
flags.transducer_receive_elevation_focus = false;

% =========================================================================
% FLAGS WHICH CAN BE CONTROLLED WITH OPTIONAL INPUTS (THESE CAN BE MODIFIED)
% =========================================================================

% flags which control the behaviour of the simulations
flags.pml_inside                = true;                 % put the PML inside the grid defined by the user
flags.record_movie              = false;                % record a movie
flags.save_to_disk              = false;                % save the input data to a HDF5 file
flags.save_to_disk_exit         = true;                 % exit the simulation after saving the HDF5 file
flags.scale_source_terms        = true;                 % apply the source scaling term to time varying sources
flags.smooth_c0                 = false;                % smooth the sound speed distribution
flags.smooth_rho0               = false;                % smooth the density distribution
flags.smooth_p0                 = true;                 % smooth the initial pressure distribution
flags.use_kspace                = true;                 % use the k-space correction
flags.use_sg                    = true;                 % use a staggered grid
flags.pml_auto                  = false;                % automatically choose the PML size to give small prime factors
flags.create_log                = false;                % create a log using diary
flags.plot_layout               = false;                % plot the simulation layout
flags.plot_sim                  = true;                 % plot the simulation
flags.plot_pml                  = true;                 % include the PML when plotting the simulation
flags.plot_scale_auto           = false;                % auto-scale the plot
flags.use_finite_difference     = false;                % use finite difference gradients instead of spectral (in 1D)
flags.mesh_plot                 = false;                % plot using mesh instead of imagesc (in 2D)
flags.stream_to_disk            = false;                % buffer the sensor data to disk (in 3D)
flags.plot_scale_log            = false;                % use a log scale for the plots
flags.data_recast               = false;                % recast the sensor data back to double precision

% =========================================================================
% VARIABLES THAT CAN BE CHANGED USING OPTIONAL INPUTS (THESE CAN BE MODIFIED)
% =========================================================================

% load the HDF5 literals (for the default compression level);
getH5Literals;

% general settings
cartesian_interp                = 'linear';              % interpolation mode for Cartesian sensor mask
hdf_compression_level           = HDF_COMPRESSION_LEVEL; % zip compression level for HDF5 input files
data_cast                       = 'off';                 % data cast
display_mask                    = 'default';             % display mask overlaid on plot
log_scale_comp_factor           = 0.02;                  % compression factor used for log scale plotting
movie_args                      = {};                    % inputs for videoWriter when saving a movie
movie_profile                   = 'Uncompressed AVI';    % compression used when saving a movie
plot_freq                       = 10;                    % number of time steps to take before updating plot
pml_search_range                = [10, 40];              % search range used when automatically determining PML size
radial_symmetry                 = 'WSWA-FFT';            % radial symmetry used in axisymmetric code
multi_axial_PML_ratio           = 0.1;                   % MPML settings

% filename for movie
if exist('MFILE', 'var')
    movie_name                  = [getDateString, '-', MFILE];
else
    movie_name                  = getDateString;
end

% default PML properties and plot scale
switch kgrid.dim
    case 1
        pml_x_alpha             = 2;
        pml_x_size              = 20;
        plot_scale              = [-1.1, 1.1];
    case 2
        pml_x_alpha             = 2;
        pml_y_alpha             = pml_x_alpha;
        pml_x_size              = 20;
        pml_y_size              = pml_x_size;
        plot_scale              = [-1, 1];
    case 3
        pml_x_alpha             = 2;
        pml_y_alpha             = pml_x_alpha;
        pml_z_alpha             = pml_x_alpha;
        pml_x_size              = 10;
        pml_y_size              = pml_x_size;
        pml_z_size              = pml_x_size;
        plot_scale              = [-1, 1];
end

% =========================================================================
% FIXED LITERALS USED IN THE CODE (THESE CAN BE MODIFIED)
% =========================================================================

% Literals used to set default parameters end with _DEF. These are cleared
% at the end of kspaceFirstOrder_inputChecking. Literals used at other
% places in the code are not cleared.

% general
STREAM_TO_DISK_STEPS_DEF        = 200;                  % number of steps before streaming to disk
COLOR_MAP                       = getColorMap;          % default color map
ESTIMATE_SIM_TIME_STEPS         = 50;                   % time steps used to estimate simulation time
HIGHEST_PRIME_FACTOR_WARNING    = 7;                    % largest prime factor before warning
KSPACE_CFL                      = 0.3;                  % default CFL value used if kgrid.t_array is set to 'auto'
PSTD_CFL                        = 0.1;                  % default CFL value used if kgrid.t_array is set to 'auto'

% source types
SOURCE_S_MODE_DEF               = 'additive';           % source mode for stress sources
SOURCE_P_MODE_DEF               = 'additive';           % source mode for pressure sources
SOURCE_U_MODE_DEF               = 'additive';           % source mode for velocity sources
DIRECTIVITY_PATTERN_DEF         = 'pressure';
DIRECTIVITY_SIZE_DEF            = 10;

% filenames
SAVE_TO_DISK_FILENAME_DEF       = 'kwave_input_data.h5';
STREAM_TO_DISK_FILENAME         = 'temp_sensor_data.bin';
LOG_NAME                        = ['k-Wave-Log-', getDateString];

% maximum scaling between p0 and plot scale before warning
switch kgrid.dim
    case 1
        PLOT_SCALE_WARNING      = 5;
    case 2
        PLOT_SCALE_WARNING      = 10;
    case 3
        PLOT_SCALE_WARNING      = 20;
end

% =========================================================================
% FIXED LITERALS USED IN THE CODE (THESE SHOULD NOT BE MODIFIED)
% =========================================================================

% number of input variables required to run the simulation codes
NUM_REQ_INPUT_VARIABLES         = 4;

% literals that link the discrete cosine and sine transform types with
% their type definitions in the functions dtt1D, dtt2D, and dtt3D
DCT1                            = 1;
DCT2                            = 2;
DCT3                            = 3;
DCT4                            = 4;
DST1                            = 5;
DST2                            = 6;
DST3                            = 7;
DST4                            = 8;