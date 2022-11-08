function test_pass = kspaceFirstOrder_1D_2D_check_checkStability(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to check the function checkStability using 1D and 2D,
%       absorbing and non-absorbing examples. 
%
% ABOUT:
%       author      - Ben Cox
%       date        - 13th August 2014
%       last update - 13th November 2014
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

% choose which tests to run
tests_to_run = [1 2 3 4];

% set up for later use
comparison1 = 0;
comparison2 = 0;
comparison3 = 0;
comparison4 = 0;

% =========================================================================
% 1D NON-ABSORBING CASE
% =========================================================================

% create the computational grid
Nx = 64;          % number of grid points in the x (row) direction
dx = 50e-3/Nx;    % grid point spacing in the x direction [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of a layered propagation medium
medium.density = 1000;          % [kg/m^3]
medium.sound_speed = 1500;      % [m/s]
medium.sound_speed_ref = 1400;  % [m/s]

% medium.sound_speed = ones(Nx,1)*1500;      % [m/s]
% medium.sound_speed(1:Nx/2,1) = 3000;
% medium.sound_speed_ref = 1400;  % [m/s]

% define a source mask
source.p_mask = zeros(Nx, 1);
source.p_mask(30) = 1;

% single source pulse
source.p = 1;

% record the final pressure field
sensor.record = {'p_final'};

% assign the input options (and
if plot_simulations
    input_args = {'PMLAlpha', 0, 'PlotSim', true , 'Smooth', false};
else
    input_args = {'PMLAlpha', 0, 'PlotSim', false, 'Smooth', false};
end

% choose the number of timepoints
Nt = 250;

if sum(tests_to_run==1)
    
    % create the time array for the stable case
    dt = 0.99*checkStability(kgrid, medium); % size of timestep [s]
    kgrid.t_array = (0:Nt)*dt;               % create array of time points

    % run the simulation for the stable case
    sensor_data_stable = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

    % create the time array for the unstable case
    dt = 1.01*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;

    % run the simulation for the unstable case
    sensor_data_unstable = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

    % the stable field should be of much lower amplitude than the unstable,
    % in which case, the test is passed. 
    comparison1 = sum(sensor_data_unstable.p_final(:).^2)/sum(sensor_data_stable.p_final(:).^2);
    if comparison1 > 1e30
        test1 = 1;
    else
        test1 = 0;
    end
    
end

% =========================================================================
% 1D ABSORBING CASE
% =========================================================================

if sum(tests_to_run==2)

    % add absorption to the medium
    medium.alpha_power = 1.5;   % [dB/(MHz^y cm)]
    medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]

    % set the reference sound speed so it does not cause the instability
    medium.sound_speed_ref = max(medium.sound_speed(:));  % [m/s]

    % create the time array and run the simulation for the stable case
    dt = 0.99*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_stable = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

    % create the time array and run the simulation for the unstable case
    dt = 1.01*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_unstable = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

    % the stable field should be of much lower amplitude than the unstable,
    % in which case the test is passed.
    comparison2 = sum(sensor_data_unstable.p_final(:).^2)/sum(sensor_data_stable.p_final(:).^2);
    if comparison2 > 1e15
        test2 = 1;
    else
        test2 = 0;
    end

end

% =========================================================================
% 2D NON-ABSORBING CASE
% =========================================================================

% create the computational grid
Nx = 64;           % number of grid points in the x (row) direction
Ny = Nx;            % number of grid points in the y (column) direction
dx = 50e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of a layered propagation medium
clear medium 
medium.density = 1000;          % [kg/m^3]
medium.sound_speed = 1500;      % [m/s]
medium.sound_speed_ref = 1400;  % [m/s]

% define a source mask
source.p_mask = zeros(Nx, Ny);
source.p_mask(30,40) = 1;

if sum(tests_to_run==3)
    
    % create the time array and run the simulation for the stable case
    dt = 0.99*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_stable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    % create the time array and run the simulation for the unstable case
    dt = 1.01*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_unstable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    % the stable field should be of much lower amplitude than the unstable,
    % in which case the test is passed.
    comparison3 = sum(sensor_data_unstable.p_final(:).^2)/sum(sensor_data_stable.p_final(:).^2);
    if comparison3 > 1e30
        test3 = 1;
    else
        test3 = 0;
    end

end

% =========================================================================
% 2D ABSORBING CASE
% =========================================================================

if sum(tests_to_run==4)

    % add absorption to the medium
    medium.alpha_power = 1.5;   % [dB/(MHz^y cm)]
    medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]

    % set the reference sound speed so it does not cause the instability
    medium.sound_speed_ref = medium.sound_speed;  % [m/s]

    % create the time array and run the simulation for the stable case
    dt = 0.99*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_stable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'UsekSpace', true);
    %sensor_data_stable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'UsekSpace', false);

    % create the time array and run the simulation for the unstable case
    dt = 1.01*checkStability(kgrid, medium);
    kgrid.t_array = (0:Nt)*dt;
    sensor_data_unstable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'UsekSpace', true);
    %sensor_data_unstable = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'UsekSpace', false);

    % the stable field should be of much lower amplitude than the unstable,
    % in which case the test is passed.
    comparison4 = sum(sensor_data_unstable.p_final(:).^2)/sum(sensor_data_stable.p_final(:).^2);
    if comparison4 > 1e15
        test4 = 1;
    else
        test4 = 0;
    end

end

% =========================================================================
% PASS THE TESTS?
% =========================================================================

% outputs the results if all the tests have been run
if (sum(tests_to_run)==10)

    if (test1+test2+test3+test4==4)
        test_pass = 1;
    else
        test_pass = 0;
    end

end

if plot_comparisons
    if comparison1>0
        display(['comparison1 = ' num2str(comparison1)])
    end
    if comparison2>0
        display(['comparison2 = ' num2str(comparison2)])
    end
    if comparison3>0
        display(['comparison3 = ' num2str(comparison3)])
    end
    if comparison4>0
        display(['comparison4 = ' num2str(comparison4)])
    end
    
end