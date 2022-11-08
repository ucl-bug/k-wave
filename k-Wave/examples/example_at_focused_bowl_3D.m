% Modelling A Focused Bowl Transducer In 3D Example
%
% This example models a focused bowl transducer in 3D. The on-axis pressure
% is compared with the exact solution calculated using focusedBowlONeil. 
%
% author: Bradley Treeby
% date: 4th February 2018
% last update: 3rd November 2022
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2022 Bradley Treeby

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

clearvars;

% =========================================================================
% DEFINE LITERALS
% =========================================================================
    
% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
model           = 1;

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 1e6;      % source frequency [Hz]
source_roc      = 30e-3;    % bowl radius of curvature [m]
source_diameter = 30e-3;    % bowl aperture diameter [m]
source_amp      = 1e6;      % source pressure [Pa]

% grid parameters
axial_size      = 50e-3;    % total grid size in the axial dimension [m]
lateral_size    = 40e-3;    % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 3;        % number of points per wavelength
t_end           = 40e-6;    % total compute time [s] (this must be long enough to reach steady state)
record_periods  = 1;        % number of periods to record
cfl             = 0.5;      % CFL number
source_x_offset = 20;       % grid points to offset the source
bli_tolerance   = 0.01;     % tolerance for truncation of the off-grid source points
upsampling_rate = 10;       % density of integration points relative to grid

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(axial_size / dx) + source_x_offset;
Ny = roundEven(lateral_size / dx);
Nz = Ny;

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period
PPP = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (PPP * source_f0);

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% --------------------
% SOURCE
% --------------------

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);

% set bowl position and orientation
bowl_pos = [kgrid.x_vec(1) + source_x_offset * kgrid.dx, 0, 0];
focus_pos = [kgrid.x_vec(end), 0, 0];

% create empty kWaveArray
karray = kWaveArray('BLITolerance', bli_tolerance, 'UpsamplingRate', upsampling_rate, 'SinglePrecision', true);

% add bowl shaped element
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);
    
% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
    
% --------------------
% MEDIUM
% --------------------

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;

% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(source_x_offset + 2:end, :, Nz/2 + 1) = 1;

% record the pressure
sensor.record = {'p'};

% record only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * PPP + 1;

% --------------------
% SIMULATION
% --------------------

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% run code
switch model
    case 1
        
        % MATLAB CPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'single', ...
            'PlotScale', [-1, 1] * source_amp);
        
    case 2
        
        % MATLAB GPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * source_amp);
        
    case 3
        
        % C++
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
        
    case 4
        
        % C++/CUDA GPU
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
        
end

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
amp = reshape(amp, Nx - source_x_offset - 1, Ny);

% extract pressure on axis
amp_on_axis = amp(:, Ny/2 + 1);

% define axis vectors for plotting
x_vec = kgrid.x_vec(source_x_offset + 2:end, :) - kgrid.x_vec(source_x_offset + 1);
y_vec = kgrid.y_vec;

% =========================================================================
% ANALYTICAL SOLUTION
% =========================================================================

% calculate the wavenumber
k = 2 * pi * source_f0 ./ c0;

% define radius and axis
x_max = (Nx * dx);
x_ref = 0:x_max/10000:x_max;

% calculate analytical solution
p_ref_axial = focusedBowlONeil(source_roc, source_diameter, source_amp / (c0 * rho0), source_f0, c0, rho0, x_ref, []);

% calculate analytical solution at exactly the same points as k-Wave
p_ref_axial_kw = focusedBowlONeil(source_roc, source_diameter, source_amp / (c0 * rho0), source_f0, c0, rho0, x_vec, []);

% calculate error
L2_error = 100 * sqrt( sum( (p_ref_axial_kw(:) - amp_on_axis(:)).^2 ) / sum( p_ref_axial_kw(:).^2 ) );
Linf_error = 100 * max(abs(p_ref_axial_kw(:) - amp_on_axis(:))) / max(p_ref_axial_kw(:));

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the pressure along the focal axis of the piston
figure;
plot(1e3 * x_ref, 1e-6 * p_ref_axial, 'k-');
hold on;
plot(1e3 * x_vec, 1e-6 * amp_on_axis, 'b.');
hold off;
set(gca, 'XLim', [0, axial_size] * 1e3);
xlabel('Axial Position [mm]');
ylabel('Pressure [MPa]');
legend('Exact', 'k-Wave', 'Location', 'Best');
title('Axial Pressure');

% plot the source mask (pml is outside the grid in this example)
figure;
subplot(1, 2, 1);
imagesc(1e3 * kgrid.y_vec, 1e3 * kgrid.x_vec, source.p_mask(:, :, ceil(Nz/2)));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
title('Source Mask');

grid_weights = karray.getArrayGridWeights(kgrid);
subplot(1, 2, 2);
imagesc(1e3 * kgrid.y_vec, 1e3 * kgrid.x_vec, grid_weights(:, :, ceil(Nz/2)));
axis image;
xlabel('y [mm]');
ylabel('x [mm]');
title('Off-Grid Source Weights');

% plot the pressure field 
figure;
imagesc(1e3 * y_vec, 1e3 * x_vec, amp);
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
axis image;
title('Pressure Field');
