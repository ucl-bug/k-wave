function test_pass = kWaveDiffusion_compare_1D_2D_3D_plane_waves(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare a plane wave in 1, 2, and 3D to catch any
%       coding bugs between the three dimensions. Eight tests are
%       performed: 
%           1. homog   + diffusion + source.T0
%           2. homog   + diffusion + source.Q
%           3. heterog + diffusion + source.T0
%           4. heterog + diffusion + source.Q
%           5. homog   + perfusion + source.T0
%           6. homog   + perfusion + source.Q
%           7. heterog + perfusion + source.T0
%           8. heterog + perfusion + source.Q
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 17th August 2015
%       last update - 7th August 2023
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2023 Bradley Treeby

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

%#ok<*UNRCH>

% set pass variable
test_pass = true;

% compare sensor data or final temperature
compare_sensor_data = false;

% set comparison thresholds
COMPARISON_THRESH_DOUBLE = 1e-13;
COMPARISON_THRESH_SINGLE = 2e-5;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 64;
Nj = 32;
dx = 1e-3;
dt = 0.5;
Nt = 100;

% define medium properties
density_1                     = 1079;     % [kg/m^3]
thermal_conductivity_1        = 0.52;     % [W/(m.K)]
specific_heat_1               = 3540;     % [J/(kg.K)]
perfusion_coeff_1             = 0.01;
blood_ambient_temperature_1   = 37;       % [degC]

density_2                     = 1200;     % [kg/m^3]
thermal_conductivity_2        = 0.13;     % [W/(m.K)]
specific_heat_2               = 4000;     % [J/(kg.K)]
perfusion_coeff_2             = 0.02;
blood_ambient_temperature_2   = 37;       % [degC]

% define source properties
ambient_temp = 37;
source_temp  = 38;
Q_amp        = 1e5;

% define source profile
source_profile = zeros(Nx, 1);
source_profile(Nx/2) = 1;
source_profile = smooth(source_profile, true);

% define optional inputs
input_args = {'PlotSim', plot_simulations, 'PlotScale', [ambient_temp, source_temp]};

% test names
test_names = {...
    'homogeneous   + diffusion + source.T0', ...
	'homogeneous   + diffusion + source.Q', ...
	'heterogeneous + diffusion + source.T0', ...
    'heterogeneous + diffusion + source.Q', ...
	'homogeneous   + perfusion + source.T0', ...
	'homogeneous   + perfusion + source.Q', ...
	'heterogeneous + perfusion + source.T0', ...
	'heterogeneous + perfusion + source.Q'};

% run tests with and without datacasting
for data_cast_ind = 2:2

    switch data_cast_ind
        case 1
            comparison_thresh = COMPARISON_THRESH_DOUBLE;
        case 2
            comparison_thresh = COMPARISON_THRESH_SINGLE;
            input_args = [input_args, {'DataCast', 'single'}]; %#ok<AGROW> 
    end

    % loop through tests
    for test_num = 1:8
    
        % clear structures
        clear source medium sensor
        
        % update command line
        disp(['Running Test: ' test_names{test_num}]);
        
        % assign medium properties 
        switch test_num
            case {1,2}
                % homog + diffusion
                medium.density                   = density_1;
                medium.thermal_conductivity      = thermal_conductivity_1;
                medium.specific_heat             = specific_heat_1;
            case {5,6}
                % homog + perfusion
                medium.density                   = density_1;
                medium.thermal_conductivity      = thermal_conductivity_1;
                medium.specific_heat             = specific_heat_1;
                medium.perfusion_coeff           = perfusion_coeff_1;
                medium.blood_ambient_temperature = blood_ambient_temperature_1;
        end
        
        % ----------------
        % 1D SIMULATION: X
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nx, dx);
    
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nx, 1);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, 1);
                medium.specific_heat             = specific_heat_1 * ones(Nx, 1);
                
                medium.density(kgrid.x > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.x > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nx, 1);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, 1);
                medium.specific_heat             = specific_heat_1 * ones(Nx, 1);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nx, 1);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nx, 1);
                
                medium.density(kgrid.x > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.x > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.x > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.x > 0) = blood_ambient_temperature_2;
        end
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * source_profile;
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * source_profile; 
        end
    
        % define sensor
        sensor.mask = ones(Nx, 1);
        
        % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_1D = kdiff.sensor_data; 
        else
            sensor_data_1D = kdiff.T;
        end
        
        % ----------------
        % 2D SIMULATION: X
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nx, dx, Nj, dx);
    
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nx, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nx, Nj);
                
                medium.density(kgrid.x > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.x > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nx, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nx, Nj);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nx, Nj);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nx, Nj);
                
                medium.density(kgrid.x > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.x > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.x > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.x > 0) = blood_ambient_temperature_2;
        end    
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * repmat(source_profile, [1, Nj]);
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * repmat(source_profile, [1, Nj]);
        end
        
        % define sensor
        sensor.mask = zeros(Nx, Nj);    
        sensor.mask(:, Nj/2) = 1;
        
        % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_2D_x = kdiff.sensor_data;
        else
            sensor_data_2D_x = reshape(kdiff.T(:, Nj/2), [Nx, 1]);
        end
        
        % ----------------
        % 2D SIMULATION: Y
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nj, dx, Nx, dx);
    
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nj, Nx);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nx);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nx);
                
                medium.density(kgrid.y > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.y > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.y > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nj, Nx);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nx);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nx);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nj, Nx);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nj, Nx);
                
                medium.density(kgrid.y > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.y > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.y > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.y > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.y > 0) = blood_ambient_temperature_2;
        end     
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * repmat(source_profile.', [Nj, 1]);
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * repmat(source_profile.', [Nj, 1]); 
        end
        
        % define sensor
        sensor.mask = zeros(Nj, Nx);    
        sensor.mask(Nj/2, :) = 1;    
        
         % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_2D_y = kdiff.sensor_data;    
        else
            sensor_data_2D_y = reshape(kdiff.T(Nj/2, :), [Nx, 1]);
        end
    
        % ----------------
        % 3D SIMULATION: X
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nx, dx, Nj, dx, Nj, dx);
        
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nx, Nj, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, Nj, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nx, Nj, Nj);
                
                medium.density(kgrid.x > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.x > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nx, Nj, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nx, Nj, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nx, Nj, Nj);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nx, Nj, Nj);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nx, Nj, Nj);
                
                medium.density(kgrid.x > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.x > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.x > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.x > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.x > 0) = blood_ambient_temperature_2;
        end
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * repmat(source_profile, [1, Nj, Nj]);
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * repmat(source_profile, [1, Nj, Nj]);
        end
        
        % define sensor
        sensor.mask = zeros(Nx, Nj, Nj);    
        sensor.mask(:, Nj/2, Nj/2) = 1;    
        
        % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_3D_x = kdiff.sensor_data;
        else
            sensor_data_3D_x = reshape(kdiff.T(:, Nj/2, Nj/2), [Nx, 1, 1]);
        end
    
        % ----------------
        % 3D SIMULATION: Y
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nj, dx, Nx, dx, Nj, dx);
    
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nj, Nx, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nx, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nx, Nj);
                
                medium.density(kgrid.y > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.y > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.y > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nj, Nx, Nj);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nx, Nj);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nx, Nj);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nj, Nx, Nj);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nj, Nx, Nj);
                
                medium.density(kgrid.y > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.y > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.y > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.y > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.y > 0) = blood_ambient_temperature_2;
        end    
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * repmat(source_profile.', [Nj, 1, Nj]);
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * repmat(source_profile.', [Nj, 1, Nj]);
        end
        
        % define sensor
        sensor.mask = zeros(Nj, Nx, Nj);    
        sensor.mask(Nj/2, :, Nj/2) = 1;    
        
        % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_3D_y = kdiff.sensor_data;
        else
            sensor_data_3D_y = reshape(kdiff.T(Nj/2, :, Nj/2), [Nx, 1, 1]);
        end
        
        % ----------------
        % 3D SIMULATION: Z
        % ----------------
    
        % create computational grid
        kgrid = kWaveGrid(Nj, dx, Nj, dx, Nx, dx);
    
        % define medium if heterog
        switch test_num
            case {3,4}
                % heterog + diffusion
                medium.density                   = density_1 * ones(Nj, Nj, Nx);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nj, Nx);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nj, Nx);
                
                medium.density(kgrid.z > 0)                 = density_2;
                medium.thermal_conductivity(kgrid.z > 0)    = thermal_conductivity_2;
                medium.specific_heat(kgrid.z > 0)           = specific_heat_2;
            case {7,8}
                % heterog + perfusion
                medium.density                   = density_1 * ones(Nj, Nj, Nx);
                medium.thermal_conductivity      = thermal_conductivity_1 * ones(Nj, Nj, Nx);
                medium.specific_heat             = specific_heat_1 * ones(Nj, Nj, Nx);
                medium.perfusion_coeff           = perfusion_coeff_1 * ones(Nj, Nj, Nx);
                medium.blood_ambient_temperature = blood_ambient_temperature_1 * ones(Nj, Nj, Nx);
                
                medium.density(kgrid.z > 0)                   = density_2;
                medium.thermal_conductivity(kgrid.z > 0)      = thermal_conductivity_2;
                medium.specific_heat(kgrid.z > 0)             = specific_heat_2;
                medium.perfusion_coeff(kgrid.z > 0)           = perfusion_coeff_2;
                medium.blood_ambient_temperature(kgrid.z > 0) = blood_ambient_temperature_2;
        end    
        
        % source
        if rem(test_num, 2)
            source.T0 = ambient_temp + (source_temp - ambient_temp) * repmat(reshape(source_profile, [1, 1, Nx]), [Nj, Nj, 1]);
        else
            source.T0 = ambient_temp;
            source.Q = Q_amp * repmat(reshape(source_profile, [1, 1, Nx]), [Nj, Nj, 1]);
        end
        
        % define sensor
        sensor.mask = zeros(Nj, Nj, Nx);    
        sensor.mask(Nj/2, Nj/2, :) = 1;    
        
        % run simulation
        kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});
        kdiff.takeTimeStep(Nt, dt);
        
        % assign output data
        if compare_sensor_data
            sensor_data_3D_z = kdiff.sensor_data;
        else
            sensor_data_3D_z = reshape(kdiff.T(Nj/2, Nj/2, :), [Nx, 1, 1]);
        end
        
        % -------------
        % COMPARISON
        % -------------
    
        diff_1D_2D_x = max(abs(sensor_data_1D(:) - sensor_data_2D_x(:)));
        diff_1D_2D_y = max(abs(sensor_data_1D(:) - sensor_data_2D_y(:)));
        diff_1D_3D_x = max(abs(sensor_data_1D(:) - sensor_data_3D_x(:)));
        diff_1D_3D_y = max(abs(sensor_data_1D(:) - sensor_data_3D_y(:)));
        diff_1D_3D_z = max(abs(sensor_data_1D(:) - sensor_data_3D_z(:)));
        
        if (diff_1D_2D_x > comparison_thresh) || ...
           (diff_1D_2D_y > comparison_thresh) || ...
           (diff_1D_3D_x > comparison_thresh) || ...
           (diff_1D_3D_y > comparison_thresh) || ...
           (diff_1D_3D_z > comparison_thresh)
            test_pass = false;
        end
        
        % -------------
        % PLOTTING
        % -------------    
        
        if plot_comparisons
            figure;
            subplot(6, 1, 1), plot(sensor_data_1D);
            title(['1D (' test_names{test_num} ')']);
            subplot(6, 1, 2), plot(sensor_data_2D_x);
            title(['2D X, L_{inf} = ' num2str(diff_1D_2D_x)]);
            subplot(6, 1, 3), plot(sensor_data_2D_y);
            title(['2D Y, L_{inf} = ' num2str(diff_1D_2D_y)]);
            subplot(6, 1, 4), plot(sensor_data_3D_x);
            title(['3D X, L_{inf} = ' num2str(diff_1D_3D_x)]);
            subplot(6, 1, 5), plot(sensor_data_3D_y);
            title(['3D Y, L_{inf} = ' num2str(diff_1D_3D_y)]);
            subplot(6, 1, 6), plot(sensor_data_3D_z);
            title(['3D Z, L_{inf} = ' num2str(diff_1D_3D_z)]);
            scaleFig(1, 1.5);
            drawnow;
        end
        
    end

end