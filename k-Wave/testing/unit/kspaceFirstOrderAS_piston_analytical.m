function test_pass = kspaceFirstOrderAS_piston_analytical(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare on-axis pressure from a circular piston source
%     simulated using k-Wave with the analytical expression from [1] using
%     the axisymmetric code.
%
%     [1] A. D. Pierce, Acoustics: An Introduction to its Physical
%     Principles and Applications. New York: Acoustical Society of America,
%     1989. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 4th February 2018
%     last update - 29th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018- Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% comparison threshold [%]
COMPARISON_THRESH = 1;

% set pass variable
test_pass = true;

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% medium parameters
C0              = 1500;     % sound speed [m/s]
RHO0            = 1000;     % density [kg/m^3]

% source parameters
SOURCE_F0       = 1e6;      % source frequency [Hz]
SOURCE_RADIUS   = 5e-3;     % piston radius [m]
SOURCE_MAG      = 1e6;      % source pressure [Pa]

% grid parameters
AXIAL_SIZE      = 32e-3;    % total grid size in the axial dimension [m]
LATERAL_SIZE    = 10e-3;    % total grid size in the lateral dimension [m]

% computational parameters
PPW             = 12;       % number of points per wavelength
T_END           = 40e-6;    % total compute time [s] (this must be long enough to reach steady state)
RECORD_PERIODS  = 3;        % number of periods to record
CFL             = 0.25;     % CFL number

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = C0 / (PPW * SOURCE_F0);   % [m]

% compute the size of the grid
Nx = roundEven(AXIAL_SIZE / dx);
Ny = roundEven(LATERAL_SIZE / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% compute points per temporal period
PPP = round(PPW / CFL);

% compute corresponding time spacing
dt = 1 / (PPP * SOURCE_F0);

% create the time array using an integer number of points per period
Nt = round(T_END / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(C0 / (dx * SOURCE_F0))]);
disp(['CFL = ' num2str(C0 * dt / dx)]);

% --------------------
% SOURCE
% --------------------

% create time varying source
source_sig = createCWSignals(kgrid.t_array, SOURCE_F0, SOURCE_MAG, 0);

% round source radius
source_radius_gp = round(SOURCE_RADIUS / dx);

% create piston source mask
source.p_mask = zeros(Nx, Ny);
source.p_mask(1, 1:source_radius_gp) = 1;

% assign source signals
source.p = source_sig;

% --------------------
% MEDIUM
% --------------------

% assign medium properties
medium.sound_speed = C0;
medium.density = RHO0;

% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny);
sensor.mask(2:end, :) = 1;

% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - RECORD_PERIODS * PPP + 1;

% --------------------
% SIMULATION
% --------------------

% run code
sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    'DataCast', 'single', ...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotSim', plot_simulations, ...
    'DisplayMask', source.p_mask, ...
    'PlotScale', [-1, 1] * SOURCE_MAG);

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, SOURCE_F0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
amp = reshape(amp, Nx - 1, Ny);

% extract pressure on axis
amp_on_axis = amp(:, 1);

% define axis vectors for plotting
x_vec = kgrid.x_vec(2:end, :) - kgrid.x_vec(1);
y_vec = kgrid.y_vec - kgrid.y_vec(1);

% =========================================================================
% ANALYTICAL SOLUTION
% =========================================================================

% calculate the wavenumber
k = 2 * pi * SOURCE_F0 ./ C0;

% define radius and axis
a = (source_radius_gp - 0.5) * dx;
x_max = (Nx * dx);
x_ref = 0:x_max/10000:x_max;

% calculate the analytical solution for a piston in an infinite baffle
% for comparison (Eq 5-7.3 in Pierce)
r = sqrt(x_ref.^2 + a^2);
p_ref = SOURCE_MAG * abs(2 * sin((k * r - k * x_ref)/2));

% get analytical solution at exactly the same points as k-Wave
r = sqrt(x_vec.^2 + a^2);
p_ref_kw = SOURCE_MAG * abs(2 * sin((k * r - k * x_vec)/2));

% calculate error
Linf_error = 100 * max(abs(p_ref_kw(:) - amp_on_axis(:))) / max(p_ref_kw(:));

% check for test pass
if Linf_error > COMPARISON_THRESH
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons

    % plot the pressure along the focal axis of the piston
    figure;
    plot(1e3 * x_ref, 1e-6 * p_ref, 'k-');
    hold on;
    plot(1e3 * x_vec, 1e-6 * amp_on_axis, 'b.');
    hold off;
    set(gca, 'XLim', [0, AXIAL_SIZE] * 1e3, 'YLim', [0, 1.1] * SOURCE_MAG * 2e-6);
    xlabel('Axial Position [mm]');
    ylabel('Pressure [MPa]');
    legend('Exact', 'k-Wave', 'Location', 'Best');
    title(['Axial Pressure (error = ' num2str(Linf_error) '%)']);

    % duplicate the pressure field
    amp = [flip(amp(:, 2:end), 2), amp];
    y_vec = [-flip(y_vec(2:end)); y_vec];

    % plot the pressure field 
    figure;
    imagesc(1e3 * x_vec, 1e3 * y_vec, amp.');
    xlabel('Axial Position [mm]');
    ylabel('Lateral Position [mm]');
    axis image;
    title('Pressure Field');
    
end