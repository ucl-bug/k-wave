% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to cast
%     loop variables to a different data type.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 26th November 2010
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2018 Bradley Treeby

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

% update command line status
disp(['  casting variables to ' data_cast ' type...']);
    
% create list of variable names used in all dimensions
if flags.elastic_code
    cast_variables = {'dt', 'mu', 'lambda'};
else
    cast_variables = {'dt', 'kappa', 'c0', 'rho0'};
end

% create a separate list for indexing variables
cast_index_variables = {};

% add variables specific to simulations in certain dimensions
switch kgrid.dim
    case 1
        
        % variables used in the fluid code
        cast_variables = [cast_variables, {...
            'ddx_k', 'shift_pos', 'shift_neg', ...
            'pml_x', 'pml_x_sgx',...
            'rho0_sgx_inv'}];
        
    case 2
        
        % variables used in both fluid and elastic codes
        cast_variables = [cast_variables, {...
            'ddx_k_shift_pos', 'ddx_k_shift_neg',...            
            'pml_x', 'pml_y',...
            'pml_x_sgx', 'pml_y_sgy',...
            'rho0_sgx_inv', 'rho0_sgy_inv'}];
        
        % y-dimension shift and derivative variables (these differ between
        % the fluid/elastic code and the axisymmetric code)
        if flags.axisymmetric
            
            % y-axis variables
            cast_variables = [cast_variables, {...
                'y_vec', 'y_vec_sg'}];
            
            % derivative and shift variables
            switch radial_symmetry
                case {'WSWA-FFT','WSWS-FFT','WS-FFT'}
                    cast_variables = [cast_variables, {...
                        'ddy_k', 'y_shift_pos', 'y_shift_neg'}];
                case 'WSWA' 
                    cast_variables = [cast_variables, {...
                        'ddy_k_wswa', 'ddy_k_hahs'}];
                case 'WSWS'
                    cast_variables = [cast_variables, {...
                        'ddy_k_wsws', 'ddy_k_haha'}];
            end

        else
            
            % derivative and shift variables in the regular code
            cast_variables = [cast_variables, {...
                'ddy_k_shift_pos', 'ddy_k_shift_neg'}];
            
        end
        
        % extra variables only used in elastic code
        if flags.elastic_code
            
            % variables used in both lossless and lossy case
            cast_variables = [cast_variables, {...
                'mu_sgxy', ...
                'mpml_x', 'mpml_y',...
                'mpml_x_sgx', 'mpml_y_sgy'}];
            
            % extra variables only used in the lossy case
            if flags.kelvin_voigt_model
                cast_variables = [cast_variables, {...
                    'chi', 'eta', 'eta_sgxy'}];                
            end
            
        end
        
    case 3
        
        % variables used in both fluid and elastic codes
        cast_variables = [cast_variables, {...
            'ddx_k_shift_pos', 'ddy_k_shift_pos', 'ddz_k_shift_pos',... 
            'ddx_k_shift_neg', 'ddy_k_shift_neg', 'ddz_k_shift_neg',...
            'pml_x', 'pml_y', 'pml_z',...
            'pml_x_sgx', 'pml_y_sgy', 'pml_z_sgz',...        
            'rho0_sgx_inv', 'rho0_sgy_inv', 'rho0_sgz_inv'}];     
        
        % extra variables only used in elastic code
        if flags.elastic_code
            
            % variables used in both lossless and lossy case
            cast_variables = [cast_variables, {...
                'mu_sgxy', 'mu_sgxz', 'mu_sgyz',...
                'mpml_x', 'mpml_y', 'mpml_z',...
                'mpml_x_sgx', 'mpml_y_sgy', 'mpml_z_sgz'}];
            
            % extra variables only used in the lossy case
            if flags.kelvin_voigt_model
                cast_variables = [cast_variables, {...
                    'chi', 'eta', 'eta_sgxy', 'eta_sgxz', 'eta_sgyz'}];                
            end
            
        end  
        
end

% add sensor mask variables
if flags.use_sensor
    cast_index_variables = [cast_index_variables, {'sensor_mask_index'}];
    if flags.binary_sensor_mask && (flags.record_u_non_staggered || flags.record_I || flags.record_I_avg)
        switch kgrid.dim
            case 1
                cast_index_variables = [cast_index_variables, {'record.x_shift_neg'}];
            case 2
                cast_index_variables = [cast_index_variables, {'record.x_shift_neg', 'record.y_shift_neg'}];
            case 3
                cast_index_variables = [cast_index_variables, {'record.x_shift_neg', 'record.y_shift_neg', 'record.z_shift_neg'}];
        end
    end
end

% additional variables only used if the medium is absorbing
if strcmp(equation_of_state, 'absorbing')
    cast_variables = [cast_variables, {'absorb_nabla1', 'absorb_nabla2', 'absorb_eta', 'absorb_tau'}];
end

% additional variables only used if the propagation is nonlinear
if flags.nonlinear
    cast_variables = [cast_variables, {'medium.BonA'}];
end

% additional variables only used if there is an initial pressure source
if flags.source_p0
    cast_variables = [cast_variables, {'source.p0'}];
end

% additional variables only used if there is a time varying pressure source
% term 
if flags.source_p
    cast_variables = [cast_variables, {'source.p'}];
    cast_index_variables = [cast_index_variables, {'p_source_pos_index'}];
    if flags.source_p_labelled
        cast_index_variables = [cast_index_variables, {'p_source_sig_index'}];
    end    
end

% additional variables only used if there is a time varying velocity source
% term 
if flags.source_ux || flags.source_uy || flags.source_uz
    cast_index_variables = [cast_index_variables, {'u_source_pos_index'}];
    if flags.source_u_labelled
        cast_index_variables = [cast_index_variables, {'u_source_sig_index'}];
    end    
end
if flags.source_ux
    cast_variables = [cast_variables, {'source.ux'}];
end
if flags.source_uy
    cast_variables = [cast_variables, {'source.uy'}];
end
if flags.source_uz
    cast_variables = [cast_variables, {'source.uz'}];
end        

% additional variables only used if there is a time varying stress source
% term 
if flags.source_sxx || flags.source_syy || flags.source_szz || flags.source_sxy || flags.source_sxz || flags.source_syz
    cast_index_variables = [cast_index_variables, {'s_source_pos_index'}];
    if flags.source_s_labelled
        cast_index_variables = [cast_index_variables, {'s_source_sig_index'}];
    end
end
if flags.source_sxx
    cast_variables = [cast_variables, {'source.sxx'}];
end
if flags.source_syy
    cast_variables = [cast_variables, {'source.syy'}];
end
if flags.source_szz
    cast_variables = [cast_variables, {'source.szz'}];
end  
if flags.source_sxy
    cast_variables = [cast_variables, {'source.sxy'}];
end
if flags.source_sxz
    cast_variables = [cast_variables, {'source.sxz'}];
end
if flags.source_syz
    cast_variables = [cast_variables, {'source.syz'}];
end

% addition variables only used if there is a transducer source
if flags.transducer_source
    cast_variables = [cast_variables, {'transducer_input_signal'}];
    cast_index_variables = [cast_index_variables, {'u_source_pos_index', 'delay_mask', 'flags.transducer_source', 'transducer_transmit_apodization'}];
end

% addition variables only used if there is a transducer sensor with an
% elevation focus
if flags.transducer_sensor && flags.transducer_receive_elevation_focus
    cast_index_variables = [cast_index_variables, {'sensor_data_buffer', 'transducer_receive_mask'}];
end

% additional variables only used with nonuniform grids
if flags.nonuniform_grid
    switch kgrid.dim
        case 1
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn'}];
        case 2
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn', 'kgrid.dyudyn'}];    
        case 3
            cast_index_variables = [cast_index_variables, {'kgrid.dxudxn', 'kgrid.dyudyn', 'kgrid.dzudzn'}];
    end
end

% additional variables only used for Cartesian sensor masks with linear
% interpolation
if flags.use_sensor && ~flags.binary_sensor_mask && ~flags.time_rev
    if  kgrid.dim == 1
        cast_variables = [cast_variables, {'record.grid_x', 'record.sensor_x'}];
    else
        cast_variables = [cast_variables, {'record.tri', 'record.bc'}];
    end
end

% additional variables only used in 2D if sensor directivity is defined
if flags.compute_directivity
    cast_variables = [cast_variables, {'sensor.directivity_angle', 'sensor.directivity_unique_angles', 'sensor.directivity_wavenumbers'}];
end

% cast variables
for cast_index = 1:length(cast_variables)
    eval([cast_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_variables{cast_index} '));']);
end

% cast index variables only if casting to the GPU
if strncmp(data_cast, 'kWaveGPU', 8)
    for cast_index = 1:length(cast_index_variables)
        eval([cast_index_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_index_variables{cast_index} '));']);
    end
end