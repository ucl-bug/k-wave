function test_pass = pstdElastic2D_check_split_field(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check that the split field components sum to give the
%     correct field, e.g., ux = ux^p + ux^s.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th November 2018
%     last update - 8th November 2018
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

% set comparison threshold
COMPARISON_THRESH   = 1e-15;

% set pass variable
test_pass = true;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% change scale to 2 to reproduce higher resolution figures in help file
scale               = 1;

% create the computational grid
PML_size            = 10;                           % [grid points]
Nx                  = 128*scale - 2*PML_size;       % [grid points]
Ny                  = 192*scale - 2*PML_size;       % [grid points]
dx                  = 0.5e-3/scale;                 % [m]
dy                  = 0.5e-3/scale;                 % [m]
kgrid               = kWaveGrid(Nx, dx, Ny, dy);

% define the medium properties for the top layer
cp1                 = 1540;   	% compressional wave speed [m/s]
cs1                 = 0;        % shear wave speed [m/s]
rho1                = 1000;     % density [kg/m^3]
alpha0_p1           = 0.1;      % compressional absorption [dB/(MHz^2 cm)]
alpha0_s1           = 0.1;      % shear absorption [dB/(MHz^2 cm)]

% define the medium properties for the bottom layer
cp2                 = 3000;     % compressional wave speed [m/s]
cs2                 = 1400;     % shear wave speed [m/s]
rho2                = 1850;     % density [kg/m^3]
alpha0_p2           = 1;        % compressional absorption [dB/(MHz^2 cm)]
alpha0_s2           = 1;        % shear absorption [dB/(MHz^2 cm)]

% create the time array
cfl                 = 0.1;
t_end               = 60e-6;
kgrid.makeTime(cp1, cfl, t_end);

% define position of heterogeneous slab
slab                = zeros(Nx, Ny);
slab(Nx/2:end, :)   = 1;

% define the source geometry in SI units (where 0, 0 is the grid center)
arc_pos             = [-15e-3, -25e-3]; % [m]
focus_pos           = [5e-3, 5e-3];     % [m]
radius              = 25e-3;            % [m]
diameter            = 20e-3;            % [m]

% define the driving signal
source_freq         = 500e3;    % [Hz]
source_strength     = 1e6;      % [Pa]
source_cycles       = 3;        % number of tone burst cycles

% define the sensor to record the maximum particle velocity everywhere
sensor.record = {'u_split_field', 'u_non_staggered'};
sensor.mask = ones(Nx, Ny);

% set the input arguments
input_args = {...
    'PlotSim', plot_simulations, ...
    'PMLSize', PML_size, ...
    'PMLAlpha', 2, ...
    'PlotPML', false, ...
    'PMLInside', false, ...
    'PlotScale', [-1, 1]*source_strength, ...
    'DisplayMask', 'off', ...
    'DataCast', 'off'};

% =========================================================================
% SIMULATION
% =========================================================================

% convert the source parameters to grid points
arc_pos     = round(arc_pos / dx) + [Nx/2, Ny/2];
focus_pos   = round(focus_pos / dx) + [Nx/2, Ny/2];
radius      = round(radius / dx);
diameter    = round(diameter / dx);

% force the diameter to be odd
if ~rem(diameter, 2)
    diameter = diameter + 1;
end

% define the medium properties
medium.sound_speed_compression            = cp1*ones(Nx, Ny);
medium.sound_speed_compression(slab == 1) = cp2;
medium.sound_speed_shear                  = cs1*ones(Nx, Ny);
medium.sound_speed_shear(slab == 1)       = cs2;
medium.density                            = rho1*ones(Nx, Ny);
medium.density(slab == 1)                 = rho2;
medium.alpha_coeff_compression            = alpha0_p1*ones(Nx, Ny);
medium.alpha_coeff_compression(slab == 1) = alpha0_p2;
medium.alpha_coeff_shear                  = alpha0_s1*ones(Nx, Ny);
medium.alpha_coeff_shear(slab == 1)       = alpha0_s2;

% generate the source geometry
source_mask = makeArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);

% assign the source
source.s_mask = source_mask;
source.sxx = -source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);
source.syy = source.sxx;

% run the elastic simulation
sensor_data_elastic = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% compute errors
diff_ux = max(abs(sensor_data_elastic.ux_non_staggered(:) - sensor_data_elastic.ux_split_p(:) - sensor_data_elastic.ux_split_s(:))) / max(abs(sensor_data_elastic.ux_non_staggered(:)));
diff_uy = max(abs(sensor_data_elastic.uy_non_staggered(:) - sensor_data_elastic.uy_split_p(:) - sensor_data_elastic.uy_split_s(:))) / max(abs(sensor_data_elastic.uy_non_staggered(:)));

% check for test pass
if (diff_ux > COMPARISON_THRESH) || (diff_uy > COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons

    % extract a single snap shot of the field
    ux   = reshape(sensor_data_elastic.ux_non_staggered(:, floor(kgrid.Nt/2)), Nx, Ny);
    uy   = reshape(sensor_data_elastic.uy_non_staggered(:, floor(kgrid.Nt/2)), Nx, Ny);
    ux_p = reshape(sensor_data_elastic.ux_split_p(:, floor(kgrid.Nt/2)), Nx, Ny);
    ux_s = reshape(sensor_data_elastic.ux_split_s(:, floor(kgrid.Nt/2)), Nx, Ny);
    uy_p = reshape(sensor_data_elastic.uy_split_p(:, floor(kgrid.Nt/2)), Nx, Ny);
    uy_s = reshape(sensor_data_elastic.uy_split_s(:, floor(kgrid.Nt/2)), Nx, Ny);

    % plot fields 
    figure;
    subplot(2, 3, 1);
    imagesc(ux);
    colorbar;
    axis image;
    title('ux');

    subplot(2, 3, 2);
    imagesc(ux_p + ux_s);
    colorbar;
    axis image;
    title('ux^p + ux^s');

    subplot(2, 3, 3);
    imagesc(abs(ux - ux_p - ux_s) ./ max(abs(ux(:))));
    colorbar;
    axis image;
    title('normalised difference');

    subplot(2, 3, 4);
    imagesc(uy);
    colorbar;
    axis image;
    title('uy');

    subplot(2, 3, 5);
    imagesc(uy_p + uy_s);
    colorbar;
    axis image;
    title('uy^p + uy^s');

    subplot(2, 3, 6);
    imagesc(abs(uy - uy_p - uy_s) ./ max(abs(uy(:))));
    colorbar;
    title('normalised difference');

    % plot components of the field
    figure;
    subplot(2, 2, 1);
    imagesc(ux_p);
    colorbar;
    axis image;
    title('ux^p');

    subplot(2, 2, 2);
    imagesc(ux_s);
    colorbar;
    axis image;
    title('ux^s');

    subplot(2, 2, 3);
    imagesc(uy_p);
    colorbar;
    axis image;
    title('uy^p');

    subplot(2, 2, 4);
    imagesc(uy_s);
    colorbar;
    axis image;
    title('uy^s');
    
end