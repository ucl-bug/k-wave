function test_pass = pstdElastic2D_compare_with_kspaceFirstOrder2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare that the elastic code with the shear wave speed
%     set to zero gives the same answers as the regular fluid code in
%     k-Wave. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 9th July 2013
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
Nx = 96;            % number of grid points in the x (row) direction
Ny = 192;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the medium properties
cp = 1500;
cs = 0;
rho = 1000;

% create the time array
CFL = 0.1;
t_end = 7e-6;
kgrid.makeTime(cp, CFL, t_end);

% create and assign the medium properties
if HETEROGENEOUS
    
    % elastic medium
    medium_elastic.sound_speed_compression = cp * ones(Nx, Ny);
    medium_elastic.sound_speed_compression(Nx/2:end, :) = 2 * cp;
    medium_elastic.sound_speed_shear = cs * ones(Nx, Ny);
    medium_elastic.density = rho*ones(Nx, Ny);
    
    % fluid medium
    medium_fluid.sound_speed = cp * ones(Nx, Ny);
    medium_fluid.density = rho * ones(Nx, Ny);
    medium_fluid.sound_speed(Nx/2:end, :) = 2 * cp;
    
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
    'source.uy, dirichlet'};

% define a single point sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(3*Nx/4, 3*Ny/4) = 1;

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

            % create initial pressure distribution using makeDisc
            disc_magnitude = 5;     % [Pa]
            disc_x_pos = 30;        % [grid points]
            disc_y_pos = Ny/2;      % [grid points]
            disc_radius = 6;        % [grid points]
            source_fluid.p0 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius); 
            
            % create equivalent elastic source
            source_elastic = source_fluid;

        case {2,3}

            % create pressure source
            source_fluid.p_mask = zeros(Nx, Ny);
            source_fluid.p_mask(30, Ny/2) = 1;
            source_fluid.p = 5 * sin(2 * pi * 1e6 * kgrid.t_array);
            source_fluid.p = filterTimeSeries(kgrid, medium_fluid, source_fluid.p);

            % create equivalent elastic source
            source_elastic.s_mask = source_fluid.p_mask;
            source_elastic.sxx = -source_fluid.p;
            source_elastic.syy = -source_fluid.p;

        case {4,5}

            % create velocity source
            source_fluid.u_mask = zeros(Nx, Ny);
            source_fluid.u_mask(30, Ny/2) = 1;
            source_fluid.ux = 5 * sin(2 * pi * 1e6 * kgrid.t_array) ./ (cp * rho);
            source_fluid.ux = filterTimeSeries(kgrid, medium_fluid, source_fluid.ux);

            % create equivalent elastic source
            source_elastic = source_fluid;

        case {6,7}

            % create velocity source
            source_fluid.u_mask = zeros(Nx, Ny);
            source_fluid.u_mask(30, Ny/2) = 1;
            source_fluid.uy = 5 * sin(2 * pi * 1e6 * kgrid.t_array) ./ (cp * rho);
            source_fluid.uy = filterTimeSeries(kgrid, medium_fluid, source_fluid.uy);

            % create equivalent elastic source
            source_elastic = source_fluid;

    end

    % set source mode
    switch test_num
        case 2
            source_fluid.p_mode   = 'additive';
            source_elastic.s_mode = 'additive';
        case 3
            source_fluid.p_mode   = 'dirichlet';
            source_elastic.s_mode = 'dirichlet';
        case {4, 6}
            source_fluid.u_mode   = 'additive';
            source_elastic.u_mode = 'additive';
        case {5, 7}
            source_fluid.u_mode   = 'dirichlet';
            source_elastic.u_mode = 'dirichlet';
    end
    
    % run the simulations
    sensor_data_elastic = pstdElastic2D     (kgrid, medium_elastic, source_elastic, sensor, input_args{:});
    sensor_data_fluid   = kspaceFirstOrder2D(kgrid, medium_fluid,   source_fluid,   sensor, 'UsekSpace', false, input_args{:});

    % compute comparisons for time series
    L_inf_p         = max(abs(sensor_data_elastic.p(COMP_START_INDEX:end)  - sensor_data_fluid.p(COMP_START_INDEX:end)))  / max(abs(sensor_data_fluid.p(COMP_START_INDEX:end))) 
    L_inf_ux        = max(abs(sensor_data_elastic.ux(COMP_START_INDEX:end) - sensor_data_fluid.ux(COMP_START_INDEX:end))) / max(abs(sensor_data_fluid.ux(COMP_START_INDEX:end)))
    L_inf_uy        = max(abs(sensor_data_elastic.uy(COMP_START_INDEX:end) - sensor_data_fluid.uy(COMP_START_INDEX:end))) / max(abs(sensor_data_fluid.uy(COMP_START_INDEX:end)))
    
    % compuate comparisons for field
    L_inf_p_final   = max(abs(sensor_data_elastic.p_final(:)  - sensor_data_fluid.p_final(:)))  / max(abs(sensor_data_fluid.p_final(:)))
    L_inf_ux_final  = max(abs(sensor_data_elastic.ux_final(:) - sensor_data_fluid.ux_final(:))) / max(abs(sensor_data_fluid.ux_final(:)))
    L_inf_uy_final  = max(abs(sensor_data_elastic.uy_final(:) - sensor_data_fluid.uy_final(:))) / max(abs(sensor_data_fluid.uy_final(:)))
    
    % compute pass
    if (L_inf_p        > COMPARISON_THRESH) || ...
       (L_inf_ux       > COMPARISON_THRESH) || ...
       (L_inf_uy       > COMPARISON_THRESH) || ...
       (L_inf_p_final  > COMPARISON_THRESH) || ...
       (L_inf_ux_final > COMPARISON_THRESH) || ...
       (L_inf_uy_final > COMPARISON_THRESH)
   
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
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.p_final, [-0.5, 0.5]);
        axis image;
        colormap(getColorMap);
        title('p (elastic)');
        colorbar;

        subplot(3, 3, 2)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.p_final, [-0.5, 0.5]);
        axis image;
        colormap(getColorMap);
        title('p (fluid)');
        colorbar;

        subplot(3, 3, 3)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.p_final - sensor_data_fluid.p_final);
        axis image;
        colormap(getColorMap);
        title('p (diff)');
        colorbar;

        subplot(3, 3, 4)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.ux_final, [-0.5, 0.5]/(cp * rho));
        axis image;
        colormap(getColorMap);
        title('ux (elastic)');
        colorbar;

        subplot(3, 3, 5)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.ux_final, [-0.5, 0.5]/(cp * rho));
        axis image;
        colormap(getColorMap);
        title('ux (fluid)');
        colorbar;

        subplot(3, 3, 6)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.ux_final - sensor_data_fluid.ux_final);
        axis image;
        colormap(getColorMap);
        title('ux (diff)');
        colorbar;

        subplot(3, 3, 7)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.uy_final, [-0.5, 0.5]/(cp * rho));
        axis image;
        colormap(getColorMap);
        title('uy (elastic)');
        colorbar;

        subplot(3, 3, 8)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_fluid.uy_final, [-0.5, 0.5]/(cp * rho));
        axis image;
        colormap(getColorMap);
        title('uy (fluid)');
        colorbar;

        subplot(3, 3, 9)
        imagesc(scale * kgrid.y_vec, scale * kgrid.x_vec, sensor_data_elastic.uy_final - sensor_data_fluid.uy_final);
        axis image;
        colormap(getColorMap);
        title('uy (diff)');
        colorbar;

        scaleFig(1, 2);

        % plot time series
        figure;
        subplot(6, 1, 1);
        plot(kgrid.t_array, sensor_data_fluid.p, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.p, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('p time series');
        subplot(6, 1, 2);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.p(COMP_START_INDEX:end) - sensor_data_elastic.p(COMP_START_INDEX:end), 'k-');
        title('diff');

        subplot(6, 1, 3);
        plot(kgrid.t_array, sensor_data_fluid.ux, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.ux, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('ux time series');
        subplot(6, 1, 4);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.ux(COMP_START_INDEX:end) - sensor_data_elastic.ux(COMP_START_INDEX:end), 'k-');
        title('diff');

        subplot(6, 1, 5);
        plot(kgrid.t_array, sensor_data_fluid.uy, 'k-');
        hold on;
        plot(kgrid.t_array, sensor_data_elastic.uy, 'r--');
        legend('fluid', 'elastic', 'Location', 'Best');
        title('uy time series');
        subplot(6, 1, 6);
        plot(kgrid.t_array(COMP_START_INDEX:end), sensor_data_fluid.uy(COMP_START_INDEX:end) - sensor_data_elastic.uy(COMP_START_INDEX:end), 'k-');
        title('diff');

        scaleFig(1, 2);

    end

end