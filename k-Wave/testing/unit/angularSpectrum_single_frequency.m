function test_pass = angularSpectrum_single_frequency(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to propagate a single frequency sine wave using
%     angularSpectrum, where the sin wave is periodically sampled, and then
%     compare the magnitude and phase of the projected signal against
%     angularSpectrumCW.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 6th March 2018
%     last update - 19th February 2019
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
COMPARISON_THRESH = 1e-12;

% set pass variable
test_pass = true;

% simulation settings
Nx = 64;
Ny = 64;
Nz = 12;
dx = 1e-4;
f0 = 2e6;
c0 = 1500;
a0 = 20;
y  = 1.5;

% calculate the time step using an integer number of points per period
cfl = 0.3;
ppw = c0 / (f0 * dx);     % points per wavelength
ppp = ceil(ppw / cfl);    % points per period
T   = 1 / f0;             % period [s]
dt  = T / ppp;            % time step [s]

% set number of points
Nt = ppp;

% create input signals
t_array = (0:Nt - 1) * dt;
sig = sin(2 * pi * f0 * t_array);

% create source shape
disc = makeDisc(Nx, Ny, Nx/2, Ny/2, Nx/4);

% create input plane
source_plane_time = repmat(reshape(sig, 1, 1, Nt), Nx, Ny, 1);
source_plane_time = bsxfun(@times, disc, source_plane_time);

% assign sound speed
medium.sound_speed = c0;

% run twice, with and without absorption
for ind = 1:2

    % if second loop, add absorption
    if ind == 2
        medium.alpha_coeff = a0;
        medium.alpha_power = y;
    end
    
    % run simulation
    [~, output_data_tv] = angularSpectrum(source_plane_time, dx, dt, (Nz - 1) * dx, medium, 'Plot', plot_simulations);

    % extract amplitude and phase
    [output_data_tv_amp, output_data_tv_phase] = extractAmpPhase( ...
        output_data_tv, 1/dt, f0, 'Dim', 3, 'FFTPadding', 1, ... 
        'Window', 'Rectangular');

    % run simulation using CW code
    output_data_cw = angularSpectrumCW(disc, dx, (Nz - 1) * dx, f0, medium);

    % extract amplitude and phase
    output_data_cw_amp = abs(output_data_cw);
    output_data_cw_phase = angle(output_data_cw);

    % account for phase offset due to retarded time
    output_data_tv_phase = output_data_tv_phase - (output_data_tv_phase(Nx/2, Nx/2) - output_data_cw_phase(Nx/2, Nx/2));
    output_data_tv_phase(output_data_tv_phase > pi) = output_data_tv_phase(output_data_tv_phase > pi) - 2 * pi;
    output_data_tv_phase(output_data_tv_phase < -pi) = output_data_tv_phase(output_data_tv_phase < -pi) + 2 * pi;

    % compute error
    err_amp = max(abs(output_data_cw_amp(:) - output_data_tv_amp(:)));
    err_phase = max(abs(output_data_cw_phase(:) - output_data_tv_phase(:)));

    % check for test pass
    if (err_amp > COMPARISON_THRESH) || (err_phase > COMPARISON_THRESH)
        test_pass = false;
    end

    % plot comparison
    if plot_comparisons

        figure;
        subplot(3, 2, 1);
        imagesc(output_data_cw_amp);
        colorbar;
        axis image
        title('CW Amp');

        subplot(3, 2, 3);
        imagesc(output_data_tv_amp);
        colorbar;
        axis image
        title('TV Amp');

        subplot(3, 2, 5);
        imagesc(abs(output_data_cw_amp - output_data_tv_amp));
        colorbar;
        axis image
        title('Difference');

        subplot(3, 2, 2);
        imagesc(output_data_cw_phase);
        colorbar;
        axis image
        title('CW Phase');

        subplot(3, 2, 4);
        imagesc(output_data_tv_phase);
        colorbar;
        axis image
        title('TV Phase');

        subplot(3, 2, 6);
        imagesc(abs(output_data_cw_phase - output_data_tv_phase));
        colorbar;
        axis image
        title('Difference');

    end
    
end