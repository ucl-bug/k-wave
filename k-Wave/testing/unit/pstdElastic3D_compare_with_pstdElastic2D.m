function test_pass = pstdElastic3D_compare_with_pstdElastic2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare an infinite line source in 2D and 3D in an
%     elastic medium to catch any coding bugs between the pstdElastic2D and
%     pstdElastic3D. 20 tests are performed:
%
%         1.  lossless + source.p0 + homogeneous
%         2.  lossless + source.p0 + heterogeneous
%         3.  lossless + source.s (additive) + homogeneous
%         4.  lossless + source.s (additive) + heterogeneous
%         5.  lossless + source.s (dirichlet) + homogeneous
%         6.  lossless + source.s (dirichlet) + heterogeneous
%         7.  lossless + source.u (additive) + homogeneous
%         8.  lossless + source.u (additive) + heterogeneous
%         9.  lossless + source.u (dirichlet) + homogeneous
%         10. lossless + source.u (dirichlet) + heterogeneous
%         11. lossy + source.p0 + homogeneous
%         12. lossy + source.p0 + heterogeneous
%         13. lossy + source.s (additive) + homogeneous
%         14. lossy + source.s (additive) + heterogeneous
%         15. lossy + source.s (dirichlet) + homogeneous
%         16. lossy + source.s (dirichlet) + heterogeneous
%         17. lossy + source.u (additive) + homogeneous
%         18. lossy + source.u (additive) + heterogeneous
%         19. lossy + source.u (dirichlet) + homogeneous
%         20. lossy + source.u (dirichlet) + heterogeneous
%
%     For each test, the infinite line source in 3D is aligned in all three
%     directions.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 13th Feb 2014
%     last update - 6th August 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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
%#ok<*NOPRT>

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set additional literals to give further permutations of the test
USE_PML             = false;
USE_SG              = true;
COMPARISON_THRESH   = 1e-14;
SMOOTH_P0_SOURCE    = false;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 64;
Ny = 64;
Nz = 32;
dx = 0.1e-3;
dy = dx;
dz = dx;

% define PML properties
PML_size    = 10;
if USE_PML
    PML_alpha   = 2;
else
    PML_alpha   = 0;
end

% define material properties
cp1      	= 1500;
cp2         = 2000;
cs1         = 0;
cs2         = 800;
rho1        = 1000;
rho2        = 1200;
alpha_p1    = 0.5;
alpha_p2    = 1;
alpha_s1    = 0.5;
alpha_s2    = 1;

% position of the heterogeneous interface
interface_position = Nx / 2;

% define time array
cfl     = 0.1;
t_end   = 3e-6;
dt      = cfl * dx / cp1;
Nt      = round(t_end / dt);
t_array = 0:dt:(Nt - 1) * dt;

% define sensor mask
sensor_mask_2D = makeCircle(Nx, Ny, Nx/2, Ny/2, 15);

% define input arguements
input_args = {'PlotScale', [-1, 1, -0.2, 0.2], 'PMLSize', PML_size, ...
    'UseSG', USE_SG, 'Smooth', false, 'PlotSim', plot_simulations};

% define source properties
source_strength     = 3;
source_position_x   = Nx/2 - 20;
source_position_y   = Ny/2 - 10;
source_freq         = 2e6;
source_signal       = source_strength * sin(2 * pi * source_freq * t_array);

% set pass variable
test_pass = true;

% test names
test_names = {...
	'lossless + source.p0 + homogeneous', ...
    'lossless + source.p0 + heterogeneous', ...
    'lossless + source.s (additive) + homogeneous', ...
    'lossless + source.s (additive) + heterogeneous', ...
    'lossless + source.s (dirichlet) + homogeneous', ...
    'lossless + source.s (dirichlet) + heterogeneous', ...
    'lossless + source.u (additive) + homogeneous', ...
    'lossless + source.u (additive) + heterogeneous', ...
    'lossless + source.u (dirichlet) + homogeneous', ...
    'lossless + source.u (dirichlet) + heterogeneous', ...
    'lossy + source.p0 + homogeneous', ...
    'lossy + source.p0 + heterogeneous', ...
    'lossy + source.s (additive) + homogeneous', ...
    'lossy + source.s (additive) + heterogeneous', ...
    'lossy + source.s (dirichlet) + homogeneous', ...
    'lossy + source.s (dirichlet) + heterogeneous', ...
    'lossy + source.u (additive) + homogeneous', ...
    'lossy + source.u (additive) + heterogeneous', ...
    'lossy + source.u (dirichlet) + homogeneous', ...
    'lossy + source.u (dirichlet) + heterogeneous', ...
    };

% lists used to set properties
p0_tests = [1, 2, 11, 12];
s_tests  = [3:6, 13:16];
u_tests  = [7:10, 17:20];
dirichlet_tests = [5, 6, 9, 10, 15, 16, 19, 20];

% set to record the velocity
sensor.record = {'u'};

% =========================================================================
% SIMULATIONS
% =========================================================================

% loop through tests
for test_num = 1:20

    % clear structures
    clear source medium

    % update command line
    disp(['Running Test: ' test_names{test_num}]);
    
    % assign medium properties
    medium.sound_speed_compression  = cp1;
    medium.sound_speed_shear        = cs1;
    medium.density                  = rho1;
    if test_num > 10
        medium.alpha_coeff_compression = alpha_p1;
        medium.alpha_coeff_shear       = alpha_s1;        
    end
        
    % ----------------
    % 2D SIMULATION
    % ----------------    
    
    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Ny, 1, 1)
    end    
    
    % source
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Ny);
        source.p0(source_position_x, source_position_y) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true); 
        end
    elseif any(s_tests == test_num)
        source.s_mask = zeros(Nx, Ny);
        source.s_mask(source_position_x, source_position_y) = 1;
        source.sxx = source_signal;
        source.syy = source_signal;
        if any(dirichlet_tests == test_num)
            source.s_mode = 'dirichlet';
        end        
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Ny);
        source.u_mask(source_position_x, source_position_y) = 1;
        source.ux = source_signal ./ (cp1 * rho1);
        source.uy = source_signal ./ (cp1 * rho1);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        end         
    else
        error('Unknown source condition.');
    end
    
    % sensor 
    sensor.mask = sensor_mask_2D;
    
    % run the simulation
    sensor_data_2D = pstdElastic2D(kgrid, medium, source, sensor, ...
        input_args{:}, 'PMLAlpha', PML_alpha);
    
    % calculate velocity amplitude
    sensor_data_2D = sqrt(sensor_data_2D.ux.^2 + sensor_data_2D.uy.^2);
    
    % ----------------
    % 3D SIMULATION: Z
    % ----------------
    
    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Ny, Nz, 1)
    end    
    
    % source
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Ny, Nz);
        source.p0(source_position_x, source_position_y, :) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true); 
        end
    elseif any(s_tests == test_num)
        source.s_mask = zeros(Nx, Ny, Nz);
        source.s_mask(source_position_x, source_position_y, :) = 1;
        source.sxx = source_signal;
        source.syy = source_signal;
        source.szz = source_signal;
        if any(dirichlet_tests == test_num)
            source.s_mode = 'dirichlet';
        end        
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Ny, Nz);
        source.u_mask(source_position_x, source_position_y, :) = 1;
        source.ux = source_signal ./ (cp1 * rho1);
        source.uy = source_signal ./ (cp1 * rho1);
        source.uz = source_signal ./ (cp1 * rho1);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        end         
    else
        error('Unknown source condition.');
    end
    
    % sensor 
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(:, :, Nz/2) = sensor_mask_2D;
    
    % run the simulation
    sensor_data_3D_z = pstdElastic3D(kgrid, medium, source, sensor, ...
        input_args{:}, 'PMLAlpha', [PML_alpha, PML_alpha, 0]);
    
    % calculate velocity amplitude
    sensor_data_3D_z = sqrt(sensor_data_3D_z.ux.^2 + sensor_data_3D_z.uy.^2);

    % ----------------
    % 3D SIMULATION: Y
    % ----------------
    
    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Nz, dz, Ny, dy);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Nz, Ny, 1)
    end    
    
    % source
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Nz, Ny);
        source.p0(source_position_x, :, source_position_y) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true); 
        end
    elseif any(s_tests == test_num)
        source.s_mask = zeros(Nx, Nz, Ny);
        source.s_mask(source_position_x, :, source_position_y) = 1;
        source.sxx = source_signal;
        source.syy = source_signal;
        source.szz = source_signal;
        if any(dirichlet_tests == test_num)
            source.s_mode = 'dirichlet';
        end        
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Nz, Ny);
        source.u_mask(source_position_x, :, source_position_y) = 1;
        source.ux = source_signal ./ (cp1 * rho1);
        source.uy = source_signal ./ (cp1 * rho1);
        source.uz = source_signal ./ (cp1 * rho1);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        end         
    else
        error('Unknown source condition.');
    end
    
    % sensor 
    sensor.mask = zeros(Nx, Nz, Ny);
    sensor.mask(:, Nz/2, :) = sensor_mask_2D;
    
    % run the simulation
    sensor_data_3D_y = pstdElastic3D(kgrid, medium, source, sensor, ...
        input_args{:}, 'PMLAlpha', [PML_alpha, 0, PML_alpha]);    
    
    % calculate velocity amplitude
    sensor_data_3D_y = sqrt(sensor_data_3D_y.ux.^2 + sensor_data_3D_y.uz.^2);    
    
    % ----------------
    % 3D SIMULATION: X
    % ----------------
    
    % create computational grid
    kgrid = kWaveGrid(Nz, dz, Nx, dx, Ny, dy);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nz, Nx, Ny, 2)
    end    
    
    % source
    if any(p0_tests == test_num)
        source.p0 = zeros(Nz, Nx, Ny);
        source.p0(:, source_position_x, source_position_y) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true); 
        end
    elseif any(s_tests == test_num)
        source.s_mask = zeros(Nz, Nx, Ny);
        source.s_mask(:, source_position_x, source_position_y) = 1;
        source.sxx = source_signal;
        source.syy = source_signal;
        source.szz = source_signal;
        if any(dirichlet_tests == test_num)
            source.s_mode = 'dirichlet';
        end        
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nz, Nx, Ny);
        source.u_mask(:, source_position_x, source_position_y) = 1;
        source.ux = source_signal ./ (cp1 * rho1);
        source.uy = source_signal ./ (cp1 * rho1);
        source.uz = source_signal ./ (cp1 * rho1);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        end         
    else
        error('Unknown source condition.');
    end
    
    % sensor 
    sensor.mask = zeros(Nz, Nx, Ny);
    sensor.mask(Nz/2, :, :) = sensor_mask_2D;
    
    % run the simulation
    sensor_data_3D_x = pstdElastic3D(kgrid, medium, source, sensor, ...
        input_args{:}, 'PMLAlpha', [0, PML_alpha, PML_alpha]);     
    
    % calculate velocity amplitude
    sensor_data_3D_x = sqrt(sensor_data_3D_x.uy.^2 + sensor_data_3D_x.uz.^2);        
    
    % -------------
    % COMPARISON
    % -------------
    
    ref_max = max(abs(sensor_data_2D(:)));
    diff_2D_3D_x = max(abs(sensor_data_2D(:) - sensor_data_3D_x(:))) / ref_max
    diff_2D_3D_y = max(abs(sensor_data_2D(:) - sensor_data_3D_y(:))) / ref_max
    diff_2D_3D_z = max(abs(sensor_data_2D(:) - sensor_data_3D_z(:))) / ref_max
    
    if (diff_2D_3D_x > COMPARISON_THRESH) || ...
       (diff_2D_3D_y > COMPARISON_THRESH) || ...
       (diff_2D_3D_z > COMPARISON_THRESH)
        test_pass = false;
    end
    
    % -------------
    % PLOTTING
    % -------------    
    
    if plot_comparisons
        
        
        figure;
        subplot(3, 3, 1)
        imagesc(sensor_data_2D);
        colorbar;
        title(['2D (' test_names{test_num} ')']);

        subplot(3, 3, 2)
        imagesc(sensor_data_3D_x);
        colorbar;
        title('3D X');

        subplot(3, 3, 3)
        imagesc(abs(sensor_data_3D_x - sensor_data_2D) ./ ref_max);
        colorbar;
        title(['L_{inf} = ' num2str(diff_2D_3D_x)]);

        subplot(3, 3, 5)
        imagesc(sensor_data_3D_y);
        colorbar;
        title('3D Y');

        subplot(3, 3, 6)
        imagesc(abs(sensor_data_3D_y - sensor_data_2D) ./ ref_max);
        colorbar;
        title(['L_{inf} = ' num2str(diff_2D_3D_y)]);
        
        subplot(3, 3, 8)
        imagesc(sensor_data_3D_z);
        colorbar;
        title('3D Z');

        subplot(3, 3, 9)
        imagesc(abs(sensor_data_3D_z - sensor_data_2D) ./ ref_max);
        colorbar;
        title(['L_{inf} = ' num2str(diff_2D_3D_z)]);
        
        figure;
        subplot(2, 1, 1);
        plot(t_array, sensor_data_2D(end/2, :));
        hold on;
        plot(t_array, sensor_data_3D_x(end/2, :));
        plot(t_array, sensor_data_3D_y(end/2, :));
        plot(t_array, sensor_data_3D_z(end/2, :));
        legend('2D', '3D x', '3D y', '3D z');
        title('Trace');
        
        subplot(2, 1, 2);
        plot(t_array, sensor_data_2D(end/2, :) - sensor_data_3D_x(end/2, :));
        hold on;
        plot(t_array, sensor_data_2D(end/2, :) - sensor_data_3D_y(end/2, :));
        plot(t_array, sensor_data_2D(end/2, :) - sensor_data_3D_z(end/2, :));
        legend('3D x', '3D y', '3D z');
        title('Error');
        
    end
    
    
end

% =========================================================================
% SUB FUNCTIONS
% =========================================================================

function setMaterialProperties(N1, N2, N3, direction)

        % sound speed and density
        medium.sound_speed_compression = cp1 * ones(N1, N2, N3);
        medium.sound_speed_shear = cs1 * ones(N1, N2, N3);
        medium.density = rho1 * ones(N1, N2, N3);
        switch direction
            case 1
                medium.sound_speed_compression(interface_position:end, :, :) = cp2;
                medium.sound_speed_shear(interface_position:end, :, :) = cs2;
                medium.density(interface_position:end, :, :) = rho2;
            case 2
                medium.sound_speed_compression(:, interface_position:end, :) = cp2;
                medium.sound_speed_shear(:, interface_position:end, :) = cs2;
                medium.density(:, interface_position:end, :) = rho2;
        end
        
        % absorption
        if isfield(medium, 'alpha_coeff_compression')
            medium.alpha_coeff_compression = alpha_p1 * ones(N1, N2, N3);
            medium.alpha_coeff_shear = alpha_s1 * ones(N1, N2, N3);
            switch direction
                case 1
                    medium.alpha_coeff_compression(interface_position:end, :, :) = alpha_p2;
                    medium.alpha_coeff_shear(interface_position:end, :, :) = alpha_s2;
                case 2
                    medium.alpha_coeff_compression(:, interface_position:end, :) = alpha_p2;
                    medium.alpha_coeff_shear(:, interface_position:end, :) = alpha_s2;
            end
            
        end    
    
end

end