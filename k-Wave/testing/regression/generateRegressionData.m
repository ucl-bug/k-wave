function generateRegressionData
% DESCRIPTION:
%     Generate data for regression tests.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 24th November 2014
%     last update - 21st October 2022
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2022 Bradley Treeby

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

% set index variable
index = 0;

% =========================================================================
% INITIAL VALUE PROBLEMS
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_ivp_homogeneous_medium';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_binary_sensor_mask';
test_examples{index}.outputs   = {'sensor_data', 'sensor_data_reordered'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_opposing_corners_sensor_mask';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_loading_external_image';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_heterogeneous_medium';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

% NOT TESTED: example_ivp_saving_movie_files (no sensor data)

index = index + 1;
test_examples{index}.name      = 'example_ivp_recording_particle_velocity';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_sensor_frequency_response';
test_examples{index}.outputs   = {'sensor_data', 'sensor_data_filtered'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_comparison_modelling_functions';
test_examples{index}.outputs   = {'sensor_data_first_order', 'sensor_data_second_order'};
test_examples{index}.precision = 'double';

% NOT TESTED: example_ivp_setting_initial_gradient (no sensor data)

index = index + 1;
test_examples{index}.name      = 'example_ivp_1D_simulation';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_3D_simulation';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_ivp_axisymmetric_simulation';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ivp_photoacoustic_waveforms';
test_examples{index}.outputs   = {'sensor_data_1D', 'sensor_data_2D', 'sensor_data_3D'};
test_examples{index}.precision = 'single';

% =========================================================================
% TIME VARYING SOURCE PROBLEMS
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_tvsp_homogeneous_medium_monopole';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_homogeneous_medium_dipole';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_transducer_field_patterns';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_steering_linear_array';
test_examples{index}.outputs   = {'source'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_snells_law';
test_examples{index}.outputs   = {'source'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_doppler_effect';
test_examples{index}.outputs   = {'sensor_data', 'approach_as', 'retreat_as'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_slit_diffraction';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_3D_simulation';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_acoustic_field_propagator';
test_examples{index}.outputs   = {'amp_out', 'phase_out'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_angular_spectrum';
test_examples{index}.outputs   = {'plane_2_kw', 'plane_2_as'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_tvsp_equivalent_source_holography';
test_examples{index}.outputs   = {'sensor_data', 'source_estimate', 'proj_asm', 'proj_dirch', 'proj_eqs'};
test_examples{index}.precision = 'single';

% =========================================================================
% SENSOR DIRECTIVITY
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_sd_focussed_detector_2D';
test_examples{index}.outputs   = {'sensor_data1', 'sensor_data2'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_sd_focussed_detector_3D';
test_examples{index}.outputs   = {'sensor_data1', 'sensor_data2'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_sd_directivity_modelling_2D';
test_examples{index}.outputs   = {'sensor_data', 'single_element_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_sd_directivity_modelling_3D';
test_examples{index}.outputs   = {'sensor_data', 'single_element_data'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_sd_sensor_directivity_2D';
test_examples{index}.outputs   = {'sensor_data1', 'sensor_data2'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_sd_directional_array_elements';
test_examples{index}.outputs   = {'sensor_data', 'element_data'};
test_examples{index}.precision = 'double';

% =========================================================================
% PHOTOACOUSTIC IMAGE RECONSTRUCTION
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_FFT_line_sensor';
test_examples{index}.outputs   = {'sensor_data', 'p_xy', 'p_xy_rs'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_pr_3D_FFT_planar_sensor';
test_examples{index}.outputs   = {'sensor_data', 'p_xyz', 'p_xyz_rs'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_TR_line_sensor';
test_examples{index}.outputs   = {'sensor_data', 'p0_recon', 'p_xy', 'p_xy_rs'};
test_examples{index}.precision = 'double';

% NOT TESTED: example_pr_2D_TR_circular_sensor (random noise added)

index = index + 1;
test_examples{index}.name      = 'example_pr_3D_TR_planar_sensor';
test_examples{index}.outputs   = {'sensor_data', 'p0_recon'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_pr_3D_TR_spherical_sensor';
test_examples{index}.outputs   = {'sensor_data', 'p0_recon', 'p0_recon_interp'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_TR_directional_sensors';
test_examples{index}.outputs   = {'sensor_data', 'p0_recon', 'p0_recon_directional'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_TR_bandlimited_sensors';
test_examples{index}.outputs   = {'sensor_data', 'p0_recon', 'p0_recon_high_pass', 'p0_recon_low_pass', 'p0_recon_gaussian'};
test_examples{index}.precision = 'double';

% NOT TESTED: example_pr_2D_TR_absorption_compensation (random noise added)
% NOT TESTED: example_pr_2D_TR_time_variant_filtering (random noise added)

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_TR_autofocus';
test_examples{index}.outputs   = {'sensor_data', 'focus_func'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_TR_iterative';
test_examples{index}.outputs   = {'sensor_data', 'p0_estimate'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_pr_2D_adjoint';
test_examples{index}.outputs   = {'sensor_data', 'p0', 'p0_1', 'p0_2', 'p0_3', 'p0_4', 'p0_5'};
test_examples{index}.precision = 'double';

% =========================================================================
% DIAGNOSTIC ULTRASOUND SIMULATION
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_us_defining_transducer';
test_examples{index}.outputs   = {'sensor_data', 'as_1', 'as_2', 'as_3'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_us_beam_patterns';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_us_transducer_as_sensor';
test_examples{index}.outputs   = {'sensor_data', 'scan_line'};
test_examples{index}.precision = 'single';

% NOT TESTED: example_us_bmode_linear_transducer (random noise added)
% NOT TESTED: example_us_bmode_phased_array (random noise added)

% =========================================================================
% NUMERICAL ANALYSIS
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_na_controlling_the_PML';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_source_smoothing';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_filtering_part_1';
test_examples{index}.outputs   = {'sensor_data', 'output_as'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_filtering_part_2';
test_examples{index}.outputs   = {'sensor_data', 'output_as'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_filtering_part_3';
test_examples{index}.outputs   = {'source_func_as', 'source_func_filtered_as', 'sensor_data', 'output_as'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_modelling_absorption';
test_examples{index}.outputs   = {'attenuation', 'attenuation_th', 'cp', 'cp_kk'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_na_modelling_nonlinearity';
test_examples{index}.outputs   = {'sensor_data', 'p_mendousse'};
test_examples{index}.precision = 'double';

% NOT TESTED: example_na_optimising_performance (no output)

% =========================================================================
% USING THE C++ CODE
% =========================================================================

% NOT TESTED

% =========================================================================
% ELASTIC WAVE PROPAGATION
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_ewp_layered_medium';
test_examples{index}.outputs   = {'sensor_data', 'sensor_data_reordered'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_ewp_plane_wave_absorption';
test_examples{index}.outputs   = {'sensor_data_comp', 'attenuation_comp', 'sensor_data_shear', 'attenuation_shear'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_ewp_shear_wave_snells_law';
test_examples{index}.outputs   = {'sensor_data_fluid', 'sensor_data_elastic'};
test_examples{index}.precision = 'single';

index = index + 1;
test_examples{index}.name      = 'example_ewp_3D_simulation';
test_examples{index}.outputs   = {'sensor_data'};
test_examples{index}.precision = 'single';

% =========================================================================
% THERMAL DIFFUSION
% =========================================================================

index = index + 1;
test_examples{index}.name      = 'example_diff_homogeneous_medium_diffusion';
test_examples{index}.outputs   = {'T_exact', 'kdiff'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_diff_homogeneous_medium_source';
test_examples{index}.outputs   = {'T_exact', 'kdiff'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_diff_binary_sensor_mask';
test_examples{index}.outputs   = {'kdiff'};
test_examples{index}.precision = 'double';

index = index + 1;
test_examples{index}.name      = 'example_diff_focused_ultrasound_heating';
test_examples{index}.outputs   = {'sensor_data', 'source', 'T1', 'T2', 'kdiff'};
test_examples{index}.precision = 'double';

% =========================================================================
% GENERATE REGRESSION DATA
% =========================================================================

% get directories
example_dir = getkWavePath('examples');
testing_dir = getkWavePath(['testing' filesep 'regression']);

% move to the example directory
cd(example_dir);

% set filename for temporary variables
temp_var_filename = [tempdir 'generate_regression_data_TEMP_VARS'];

% clear previous regression test results (this name is used explicitly
% later in the file)
if exist(temp_var_filename, 'file')
    delete(temp_var_filename);
end

% run the examples one by one
for filename_index = 1:length(test_examples)
    
    % move to the example directory
    cd(example_dir);
    
    % save required variables to the workspace to circumvent clear all
    % within the example files
    save([tempdir 'generate_regression_data_TEMP_VARS'], 'example_dir', 'testing_dir', 'filename_index', 'test_examples');
    
    % get the current m-file name
    fn = test_examples{filename_index}.name;
    
    % display the filename
    disp(' ');
    disp(['Saving regression data for ' fn ' (Test ' num2str(filename_index) ' of ' num2str(length(test_examples)) ')']);

    % run the example file
    run([example_dir fn])
    
    % close the figure windows
    close all;
    
    % load back the workspace variables cleared by the example file
    load([tempdir 'generate_regression_data_TEMP_VARS']);
    
    % move to the testing directory
    cd(testing_dir);
    
    % generate computer info and add precision
    comp_info = getComputerInfo;
    comp_info.precision = test_examples{filename_index}.precision;
    
    % save the regression data
    save(test_examples{filename_index}.name, test_examples{filename_index}.outputs{:}, 'comp_info')
    
end

% remove temp data
delete([tempdir 'generate_regression_data_TEMP_VARS']);