function test_pass = pstdElastic3D_check_mpml_stability(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to test the stability of the pml and m-pml.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 12 August 2014
%       last update - 12 August 2014
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

% set pass variable
test_pass = true;

% create the computational grid
PML_SIZE = 10;
Nx = 80 - 2*PML_SIZE;
Ny = 64 - 2*PML_SIZE;
Nz = 64 - 2*PML_SIZE;
dx = 0.1e-3;          
dy = 0.1e-3;          
dz = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500*ones(Nx, Ny, Nz); % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny, Nz);     % [m/s]
medium.density                 = 1000*ones(Nx, Ny, Nz); % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :, :) = 2000; % [m/s]
medium.sound_speed_shear(Nx/2:end, :, :)       = 1000;  % [m/s]
medium.density(Nx/2:end, :, :)                 = 1200; % [kg/m^3]

% create the time array
cfl   = 0.3;        % Courant-Friedrichs-Lewy number
t_end = 8e-6;      % [s]
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed_compression(:)), cfl, t_end);

% define the source mask
s_rad = 15;
s_height = 8;
offset = 15;
ss = makeSphericalSection(s_rad, s_height);
ss_width = size(ss, 2);
ss_half_width = floor(ss_width/2);
y_start_pos = Ny/2 - ss_half_width;
y_end_pos = y_start_pos + ss_width - 1;
z_start_pos = Nz/2 - ss_half_width;
z_end_pos = z_start_pos + ss_width - 1;
source.s_mask = zeros(Nx, Ny, Nz);
source.s_mask(1+offset:s_height+offset, y_start_pos:y_end_pos, z_start_pos:z_end_pos) = ss;
source.s_mask(:, :, Nz/2:end) = 0;

% define the source signal
source.sxx = toneBurst(1/kgrid.dt, 1e6, 3);
source.syy = source.sxx;
source.szz = source.sxx;

% define sensor
sensor.record = {'u_final'};
sensor.mask = [];

% define input arguments
input_args = {'PlotScale', [-0.1, 0.1, -0.005, 0.005], 'DataCast', 'single',...
    'DisplayMask', 'off', 'PMLAlpha', 2, 'PMLInside', false,...
    'PMLSize', PML_SIZE, 'PlotSim', plot_simulations};

% run the simulations
sensor_data_pml = pstdElastic3D(kgrid, medium, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0);
sensor_data_mpml = pstdElastic3D(kgrid, medium, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0.1);

% check magnitudes
pml_max = max([sensor_data_pml.ux_final(:); sensor_data_pml.uy_final(:); sensor_data_pml.uz_final(:)]);
mpml_max = max([sensor_data_mpml.ux_final(:); sensor_data_mpml.uy_final(:); sensor_data_mpml.uz_final(:)]);

% set reference magnitude (initial source)
ref_max = 1 / ( max(medium.sound_speed_shear(:)) * max(medium.density(:)) );

% check results - the test should fail if the pml DOES work (i.e., it
% doesn't become unstable), or if the m-pml DOESN'T work (i.e., it does
% become unstable)
if pml_max < ref_max || mpml_max > ref_max
    test_pass = false;
end

% plot if required
if plot_comparisons
    figure;
    subplot(2, 3, 1)
    imagesc(squeeze(sensor_data_pml.ux_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_x with PML');
    
    subplot(2, 3, 2)
    imagesc(squeeze(sensor_data_pml.uy_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_y with PML');
    
    subplot(2, 3, 3)
    imagesc(squeeze(sensor_data_pml.uz_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_z with PML');    
    
    subplot(2, 3, 4)
    imagesc(squeeze(sensor_data_mpml.ux_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_x with M-PML');
    
    subplot(2, 3, 5)
    imagesc(squeeze(sensor_data_mpml.uy_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_y with M-PML');
    
    subplot(2, 3, 6)
    imagesc(squeeze(sensor_data_mpml.uz_final(:, :, end/2)));
    axis image;
    colorbar;
    title('u_z with M-PML');    
end