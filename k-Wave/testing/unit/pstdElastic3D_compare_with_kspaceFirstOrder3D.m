function test_pass = pstdElastic3D_compare_with_kspaceFirstOrder3D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare that the elastic code with the shear wave speed
%     set to zero gives the same answers as the regular fluid code in
%     k-Wave. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 16th July 2013
%     last update - 18th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2017 Bradley Treeby

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
HETEROGENEOUS       = true;
USE_PML             = false;
DATA_CAST           = 'off';
COMPARISON_THRESH   = 5e-13;

% option to skip the first point in the time series (for p0 sources, there
% is a strange bug where there is a high error for the first stored time
% point)
COMP_START_INDEX    = 2;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64; 
Ny = 64; 
Nz = 64; 
dx = 0.1e-3;
dy = 0.1e-3;
dz = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the medium properties
cp = 1500;
cs = 0;
rho = 1000;

% create the time zarray
CFL = 0.1;
t_end = 3e-6;
kgrid.makeTime(cp, CFL, t_end);

% create and assign the variables
if HETEROGENEOUS
    
    % elastic medium
    medium_elastic.sound_speed_compression = cp * ones(Nx, Ny, Nz);
    medium_elastic.sound_speed_shear = cs * ones(Nx, Ny, Nz);
    medium_elastic.density = rho * ones(Nx, Ny, Nz);
    medium_elastic.sound_speed_compression(Nx/2:end, :, :) = 2 * cp;
    
    % fluid medium
    medium_fluid.sound_speed = cp * ones(Nx, Ny, Nz);
    medium_fluid.density = rho * ones(Nx, Ny, Nz);
    medium_fluid.sound_speed(Nx/2:end, :, :) = 2 * cp;
    
else
    
    % elastic medium
    medium_elastic.sound_speed_compression = cp;
    medium_elastic.sound_speed_shear = cs;
    medium_elastic.density = rho;
    
    % fluid medium
    medium_fluid.sound_speed = cp;
    medium_fluid.density = rho;
    
end

% set pass variable
test_pass = true;

% test names
test_names = {...
    'source.p0', ...
    'source.p, additive', ...
    'source.p, dirichlet', ...
    'source.ux, additive', ...
    'source.ux, dirichlet', ...
    'source.uy, additive', ...
    'source.uy, dirichlet',...
    'source.uz, additive', ...
    'source.uz, dirichlet'};

% define a single point sensor
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(3*Nx/4, 3*Ny/4, 3*Nz/4) = 1;

% set some things to record
sensor.record = {'p', 'p_final', 'u', 'u_final'};

% set input args
input_args = {'PlotSim', plot_simulations, 'DataCast', DATA_CAST};

if ~USE_PML
    input_args = [input_args {'PMLAlpha', 0}]; 
end

% loop through tests
for test_num = 1:7

    % clear structures
    clear source_fluid source_elastic
    
    % update command line
    disp(['Running Test: ' test_names{test_num}]);
    
    switch test_num

        case 1

            % create initial pressure distribution using makeBall
            disc_magnitude = 5;         % [Pa]
            disc_x_pos = Nx/2 - 10;     % [grid points]
            disc_y_pos = Ny/2;          % [grid points]
            disx_z_pos = Nz/2;          % [grid points]
            disc_radius = 3;            % [grid points]
            source_fluid.p0 = disc_magnitude * makeBall(Nx, Ny, Nz, disc_x_pos, disc_y_pos, disx_z_pos, disc_radius);

            % assign to elastic source
            source_elastic = source_fluid;

        case {2,3}

            % create pressure source
            source_fluid.p_mask = zeros(Nx, Ny, Nz);
            source_fluid.p_mask(Nx/2 - 10, Ny/2, Nz/2) = 1;
            source_fluid.p = 20 * sin(2 * pi * 1e6 * kgrid.t_array);
            source_fluid.p = filterTimeSeries(kgrid, medium_fluid, source_fluid.p);

            % create equivalent stress source
            source_elastic.s_mask = zeros(Nx, Ny, Nz);
            source_elastic.s_mask(Nx/2 - 10, Ny/2, Nz/2) = 1;
            source_elastic.sxx = -source_fluid.p;
            source_elastic.syy = -source_fluid.p;
            source_elastic.szz = -source_fluid.p;
      
        case {4,5}

            % create velocity source
            source_fluid.u_mask = zeros(Nx, Ny, Nz);
            source_fluid.u_mask(Nx/2 - 10, Ny/2, Nz/2) = 1;
            source_fluid.ux = 20 * sin(2 * pi * 1e6 * kgrid.t_array) ./ (cp * rho);
            source_fluid.ux = filterTimeSeries(kgrid, medium_fluid, source_fluid.ux);

            % assign to elastic source
            source_elastic = source_fluid;

        case {6,7}

            % create velocity source
            source_fluid.u_mask = zeros(Nx, Ny, Nz);
            source_fluid.u_mask(Nx/2 - 10, Ny/2, Nz/2) = 1;
            source_fluid.uy = 20 * sin(2 * pi * 1e6 * kgrid.t_array) ./ (cp * rho);
            source_fluid.uy = filterTimeSeries(kgrid, medium_fluid, source_fluid.uy);

            % assign to elastic source
            source_elastic = source_fluid;

        case {8,9}

            % create velocity source
            source_fluid.u_mask = zeros(Nx, Ny, Nz);
            source_fluid.u_mask(Nx/2 - 10, Ny/2, Nz/2) = 1;
            source_fluid.uz = 20 * sin(2 * pi * 1e6 * kgrid.t_array) ./ (cp * rho);
            source_fluid.uz = filterTimeSeries(kgrid, medium_fluid, source_fluid.uz);

            % assign to elastic source
            source_elastic = source_fluid;
        
    end

    % set source mode
    switch test_num
        case 2
            source_fluid.p_mode = 'additive';
            source_elastic.s_mode = 'additive';
        case 3
            source_fluid.p_mode = 'dirichlet';
            source_elastic.s_mode = 'dirichlet';
        case {4, 6}
            source_fluid.u_mode = 'additive';
            source_elastic.u_mode = 'additive';
        case {5, 7}
            source_fluid.u_mode = 'dirichlet';
            source_elastic.u_mode = 'dirichlet';
    end

    % run the simulations
    sensor_data_elastic = pstdElastic3D(kgrid, medium_elastic, source_elastic, sensor, input_args{:});
    sensor_data_fluid = kspaceFirstOrder3D(kgrid, medium_fluid, source_fluid, sensor, input_args{:}, 'UsekSpace', false);

    % compute comparisons for time series
    L_inf_p        = max(abs(sensor_data_elastic.p(COMP_START_INDEX:end)  - sensor_data_fluid.p(COMP_START_INDEX:end)))  / max(abs(sensor_data_fluid.p(COMP_START_INDEX:end)))
    L_inf_ux       = max(abs(sensor_data_elastic.ux(COMP_START_INDEX:end) - sensor_data_fluid.ux(COMP_START_INDEX:end))) / max(abs(sensor_data_fluid.ux(COMP_START_INDEX:end)))
    L_inf_uy       = max(abs(sensor_data_elastic.uy(COMP_START_INDEX:end) - sensor_data_fluid.uy(COMP_START_INDEX:end))) / max(abs(sensor_data_fluid.uy(COMP_START_INDEX:end)))
    L_inf_uz       = max(abs(sensor_data_elastic.uz(COMP_START_INDEX:end) - sensor_data_fluid.uz(COMP_START_INDEX:end))) / max(abs(sensor_data_fluid.uz(COMP_START_INDEX:end)))
    
    % compuate comparisons for field
    L_inf_p_final  = max(abs(sensor_data_elastic.p_final(:)  - sensor_data_fluid.p_final(:)))  / max(abs(sensor_data_fluid.p_final(:)))
    L_inf_ux_final = max(abs(sensor_data_elastic.ux_final(:) - sensor_data_fluid.ux_final(:))) / max(abs(sensor_data_fluid.ux_final(:)))
    L_inf_uy_final = max(abs(sensor_data_elastic.uy_final(:) - sensor_data_fluid.uy_final(:))) / max(abs(sensor_data_fluid.uy_final(:)))
    L_inf_uz_final = max(abs(sensor_data_elastic.uz_final(:) - sensor_data_fluid.uz_final(:))) / max(abs(sensor_data_fluid.uz_final(:)))
    
    % compute pass
    if (L_inf_p        > COMPARISON_THRESH) || ...
       (L_inf_ux       > COMPARISON_THRESH) || ...
       (L_inf_uy       > COMPARISON_THRESH) || ...
       (L_inf_uz       > COMPARISON_THRESH) || ...
       (L_inf_p_final  > COMPARISON_THRESH) || ...
       (L_inf_ux_final > COMPARISON_THRESH) || ...
       (L_inf_uy_final > COMPARISON_THRESH) || ...
       (L_inf_uz_final > COMPARISON_THRESH)

        % set test variable
        test_pass = false;
        
        % display result
        disp('Test Component Failed!');   
        
    else
        
        % display result
        disp('Test Component Passed.');        
   
    end
    
    % =========================================================================
    % PLOT COMPARISONS
    % =========================================================================

    if plot_comparisons

        [x_sc, scale] = scaleSI(kgrid.x_vec); %#ok<ASGLU>

        % plot final fields
        figure;
        subplot(3, 3, 1)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.p_final(:, :, end/2), [-0.5, 0.5]);
        axis image;
        colormap(getColorMap);
        title('p (elastic)');
        colorbar;

        subplot(3, 3, 2)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.p_final(:, :, end/2), [-0.5, 0.5]);
        axis image;
        colormap(getColorMap);
        title('p (fluid)');
        colorbar;

        subplot(3, 3, 3)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.p_final(:, :, end/2) - sensor_data_fluid.p_final(:, :, end/2));
        axis image;
        colormap(getColorMap);
        title('p (diff)');
        colorbar;

        subplot(3, 3, 4)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.ux_final(:, :, end/2), [-0.5, 0.5]/(cp*rho));
        axis image;
        colormap(getColorMap);
        title('ux (elastic)');
        colorbar;

        subplot(3, 3, 5)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.ux_final(:, :, end/2), [-0.5, 0.5]/(cp*rho));
        axis image;
        colormap(getColorMap);
        title('ux (fluid)');
        colorbar;

        subplot(3, 3, 6)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.ux_final(:, :, end/2) - sensor_data_fluid.ux_final(:, :, end/2));
        axis image;
        colormap(getColorMap);
        title('ux (diff)');
        colorbar;

        subplot(3, 3, 7)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.uy_final(:, :, end/2), [-0.5, 0.5]/(cp*rho));
        axis image;
        colormap(getColorMap);
        title('uy (elastic)');
        colorbar;

        subplot(3, 3, 8)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.uy_final(:, :, end/2), [-0.5, 0.5]/(cp*rho));
        axis image;
        colormap(getColorMap);
        title('uy (fluid)');
        colorbar;

        subplot(3, 3, 9)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.uy_final(:, :, end/2) - sensor_data_fluid.uy_final(:, :, end/2));
        axis image;
        colormap(getColorMap);
        title('uy (diff)');
        colorbar;

        scaleFig(1, 2);

        % plot time series
        h = figure;
        subplot(4, 2, 1);
        plot(kgrid.t_array, sensor_data_fluid.p, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.p, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('p time series');
        subplot(4, 2, 2);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.p(COMP_START_INDEX:end) - sensor_data_elastic.p(COMP_START_INDEX:end), 'k-');
        title('diff');

        subplot(4, 2, 3);
        plot(kgrid.t_array, sensor_data_fluid.ux, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.ux, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('ux time series');
        subplot(4, 2, 4);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.ux(COMP_START_INDEX:end) - sensor_data_elastic.ux(COMP_START_INDEX:end), 'k-');
        title('diff');

        subplot(4, 2, 5);
        plot(kgrid.t_array, sensor_data_fluid.uy, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.uy, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('uy time series');
        subplot(4, 2, 6);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.uy(COMP_START_INDEX:end) - sensor_data_elastic.uy(COMP_START_INDEX:end), 'k-');
        title('diff');

        subplot(4, 2, 7);
        plot(kgrid.t_array, sensor_data_fluid.uz, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.uz, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('uz time series');
        subplot(4, 2, 8);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.uz(COMP_START_INDEX:end) - sensor_data_elastic.uz(COMP_START_INDEX:end), 'k-');
        title('diff');

        scaleFig(2, 2);
        drawnow;
        refresh(h);
        
    end
end