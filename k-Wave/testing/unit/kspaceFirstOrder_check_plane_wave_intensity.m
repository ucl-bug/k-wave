function test_pass = kspaceFirstOrder_check_plane_wave_intensity(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to check the intensity of a plane wave in 1, 2, and 3D
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 14th July 2014
%       last update - 14th July 2014
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

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

%#ok<*NOPRT>
%#ok<*UNRCH>

% set additional literals to give further permutations of the test
STAGGERED_GRID      = true;
USE_KSPACE          = true;
COMPARISON_THRESH   = 1;        % [percentage]

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 64;
Nj = 32;
dx = 0.05e-3;

% define PML properties
PML_size = 20;
PML_alpha = 2;
    
% define medium properties
medium.sound_speed = 1500;
medium.density = 1000;

% set the source parameters
source_freq = 10e6;      % [Hz]
source_intensity = 10;   % [W/cm^2] time averaged intensity
source_strength = sqrt(2*source_intensity*1e4*medium.density*medium.sound_speed);

% set which variables to record
sensor.record = {'I', 'I_avg'};

% calculate the time length of one acoustic period
T = 1/source_freq;

% create the time array using a nice number of points
points_per_period = 30;
dt = T/(2*points_per_period);
t_end = 3e-6;
t_array = 0:dt:t_end;
Nt = length(t_array);

% set start point
sensor.record_start_index = Nt - 2*points_per_period + 1;

% create the source
source.p = source_strength*sin(2*pi*source_freq*t_array);

% smooth startup using a cosine ramp
ramp_length = 200;
ramp = (-cos( (0:(ramp_length-1))*pi/ramp_length ) + 1)/2;
source.p(1:length(ramp)) = source.p(1:length(ramp)) .* ramp;

% double check the CFL
disp(['CFL = ' num2str(dt*medium.sound_speed/dx)]);

% define source and sensor properties
source_position = PML_size + 1;
sensor_position = Nx - PML_size - 1; 

% set pass variable
test_pass = true;

% ----------------
% 1D SIMULATION: X
% ----------------

% create computational grid
kgrid = kWaveGrid(Nx, dx);
kgrid.t_array = t_array;

% define optional inputs
input_args = {'PlotScale', 'auto', ...
    'UsekSpace', USE_KSPACE, 'UseSG', STAGGERED_GRID, ...
    'PlotSim', plot_simulations};

% source
source.p_mask = zeros(Nx, 1);
source.p_mask(source_position) = 1;

% sensor
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_position) = 1;

% run simulation
sensor_data_1D_x = kspaceFirstOrder1D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size);

% ----------------
% 2D SIMULATION: X
% ----------------

% create computational grid
kgrid = kWaveGrid(Nx, dx, Nj, dx);
kgrid.t_array = t_array;

% source
source.p_mask = zeros(Nx, Nj);
source.p_mask(source_position, :) = 1; 

% sensor
sensor.mask = zeros(Nx, Nj);
sensor.mask(sensor_position, Nj/2) = 1;

% run simulation
sensor_data_2D_x = kspaceFirstOrder2D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', [PML_alpha 0], 'PMLSize', [PML_size 2]);

% ----------------
% 2D SIMULATION: Y
% ----------------

% create computational grid
kgrid = kWaveGrid(Nj, dx, Nx, dx);
kgrid.t_array = t_array;

% source
source.p_mask = zeros(Nj, Nx);
source.p_mask(:, source_position) = 1; 

% sensor
sensor.mask = zeros(Nj, Nx);
sensor.mask(Nj/2, sensor_position) = 1;

% run simulation
sensor_data_2D_y = kspaceFirstOrder2D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', [0 PML_alpha], 'PMLSize', [2 PML_size]);    

% ----------------
% 3D SIMULATION: X
% ----------------

% create computational grid
kgrid = kWaveGrid(Nx, dx, Nj, dx, Nj, dx);
kgrid.t_array = t_array;

% source
source.p_mask = zeros(Nx, Nj, Nj);
source.p_mask(source_position, :, :) = 1;      

% sensor
sensor.mask = zeros(Nx, Nj, Nj);
sensor.mask(sensor_position, Nj/2, Nj/2) = 1;

% run simulation
sensor_data_3D_x = kspaceFirstOrder3D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', [PML_alpha 0 0], 'PMLSize', [PML_size 2 2]);

 % ----------------
% 3D SIMULATION: Y
% ----------------

% create computational grid
kgrid = kWaveGrid(Nj, dx, Nx, dx, Nj, dx);
kgrid.t_array = t_array;

% source
source.p_mask = zeros(Nj, Nx, Nj);
source.p_mask(:, source_position, :) = 1;      

% sensor
sensor.mask = zeros(Nj, Nx, Nj);
sensor.mask(Nj/2, sensor_position, Nj/2) = 1;

% run simulation
sensor_data_3D_y = kspaceFirstOrder3D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', [0 PML_alpha 0], 'PMLSize', [2 PML_size 2]);

 % ----------------
% 3D SIMULATION: Z
% ----------------

% create computational grid
kgrid = kWaveGrid(Nj, dx, Nj, dx, Nx, dx);
kgrid.t_array = t_array;

% source
source.p_mask = zeros(Nj, Nj, Nx);
source.p_mask(:, :, source_position) = 1;      

% sensor
sensor.mask = zeros(Nj, Nj, Nx);
sensor.mask(Nj/2, Nj/2, sensor_position) = 1;

% run simulation
sensor_data_3D_z = kspaceFirstOrder3D(kgrid, medium, source, sensor,...
    input_args{:}, 'PMLAlpha', [0 0 PML_alpha], 'PMLSize', [2 2 PML_size]);

% -------------
% COMPARISON
% -------------

% convert k-Wave intensity from W/m^2 to W/cm^2 and compare with desired
% intensity
error_intensity_1D_x = 100*abs(source_intensity - 1e-4*sensor_data_1D_x.Ix_avg)/source_intensity;
error_intensity_2D_x = 100*abs(source_intensity - 1e-4*sensor_data_2D_x.Ix_avg)/source_intensity;
error_intensity_2D_y = 100*abs(source_intensity - 1e-4*sensor_data_2D_y.Iy_avg)/source_intensity;
error_intensity_3D_x = 100*abs(source_intensity - 1e-4*sensor_data_3D_x.Ix_avg)/source_intensity;
error_intensity_3D_y = 100*abs(source_intensity - 1e-4*sensor_data_3D_y.Iy_avg)/source_intensity;
error_intensity_3D_z = 100*abs(source_intensity - 1e-4*sensor_data_3D_z.Iz_avg)/source_intensity;

if (error_intensity_1D_x > COMPARISON_THRESH) || ...
   (error_intensity_2D_x > COMPARISON_THRESH) || ...
   (error_intensity_2D_y > COMPARISON_THRESH) || ...
   (error_intensity_3D_x > COMPARISON_THRESH) || ...
   (error_intensity_3D_y > COMPARISON_THRESH) || ...
   (error_intensity_3D_z > COMPARISON_THRESH)
    test_pass = false;
end

% -------------
% PLOTTING
% -------------    

if plot_comparisons
    figure;
    subplot(6, 1, 1), plot(sensor_data_1D_x.Ix);
    title(['1D X, I_{avg} = ' num2str(sensor_data_1D_x.Ix_avg*1e-4) ' W/cm^2']);
    subplot(6, 1, 2), plot(sensor_data_2D_x.Ix);
    title(['2D X, I_{avg} = ' num2str(sensor_data_2D_x.Ix_avg*1e-4) ' W/cm^2']);
    subplot(6, 1, 3), plot(sensor_data_2D_y.Iy);
    title(['2D Y, I_{avg} = ' num2str(sensor_data_2D_y.Iy_avg*1e-4) ' W/cm^2']);
    subplot(6, 1, 4), plot(sensor_data_3D_x.Ix);
    title(['3D X, I_{avg} = ' num2str(sensor_data_3D_x.Ix_avg*1e-4) ' W/cm^2']);
    subplot(6, 1, 5), plot(sensor_data_3D_y.Iy);
    title(['3D Y, I_{avg} = ' num2str(sensor_data_3D_y.Iy_avg*1e-4) ' W/cm^2']);
    subplot(6, 1, 6), plot(sensor_data_3D_z.Iz);
    title(['3D Z, I_{avg} = ' num2str(sensor_data_3D_z.Iz_avg*1e-4) ' W/cm^2']);
    scaleFig(1, 1.5);
    drawnow;
end