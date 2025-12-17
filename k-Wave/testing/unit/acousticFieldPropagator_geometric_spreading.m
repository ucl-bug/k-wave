function test_pass = acousticFieldPropagator_geometric_spreading(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to compare geometric spreading with analytical values.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 9th January 2019
%     last update - 9th January 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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
end

% set comparison threshold
comparison_threshold = 0.01;

% set pass variable
test_pass = true;

% =========================================================================
% SETUP
% =========================================================================

% literals
USE_RAMP = true;
TIME_EXP_FAC = 1.1;

% define grid size
Nx = 200;
Ny = 64;
Nz = 64;
dx = 0.1e-3;

% define source frequency and sound speed
c0 = 1500;
f = 1e6;
input_phase = 0;

% create distance axis
x = (0:Nx - 1) * dx;

% normalise to 1 at a distance of 1 mm
[~, ind] = findClosest(x, 1e-3);

% =========================================================================
% 1D
% =========================================================================

% define point source
aperture_amp = zeros(Nx, 1);
aperture_amp(1) = 1;

% calculate pressure
[amp, ~] = acousticFieldPropagator(aperture_amp, input_phase, dx, f, c0, 'UseRamp', USE_RAMP, 'TimeExpansionFactor', TIME_EXP_FAC);

% extract the radial decay
amp_dec_1D = amp;

% normalise
amp_dec_1D = amp_dec_1D ./ amp_dec_1D(ind);

% calculate analytical decay
amp_dec_1D_an = ones(Nx, 1);

% compute error
err = max(abs(amp_dec_1D(ind:end) - amp_dec_1D_an(ind:end)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end
    
% =========================================================================
% 2D
% =========================================================================

% define point source
aperture_amp = zeros(Nx, Ny);
aperture_amp(1, Ny/2) = 1;

% calculate pressure
[amp, ~] = acousticFieldPropagator(aperture_amp, input_phase, dx, f, c0, 'UseRamp', USE_RAMP, 'TimeExpansionFactor', TIME_EXP_FAC);

% extract the radial decay
amp_dec_2D = amp(1:end, Ny/2);

% calculate analytical decay
amp_dec_2D_an = (1 ./ sqrt(x)).';

% normalise
amp_dec_2D = amp_dec_2D ./ amp_dec_2D(ind);
amp_dec_2D_an = amp_dec_2D_an ./ amp_dec_2D_an(ind);

% compute error
err = max(abs(amp_dec_2D(ind:end) - amp_dec_2D_an(ind:end)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end
    
% =========================================================================
% 3D
% =========================================================================

% define point source
aperture_amp = zeros(Nx, Ny, Nz);
aperture_amp(1, Ny/2, Nz/2) = 1;

% calculate pressure
[amp, ~] = acousticFieldPropagator(aperture_amp, input_phase, dx, f, c0, 'UseRamp', USE_RAMP, 'TimeExpansionFactor', TIME_EXP_FAC);

% extract the radial decay
amp_dec_3D = amp(1:end, Ny/2, Nz/2);

% calculate analytical decay
amp_dec_3D_an = (1 ./ x).';

% normalise
amp_dec_3D = amp_dec_3D ./ amp_dec_3D(ind);
amp_dec_3D_an = amp_dec_3D_an ./ amp_dec_3D_an(ind);

% compute error
err = max(abs(amp_dec_3D(ind:end) - amp_dec_3D_an(ind:end)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end
    
% =========================================================================
% PLOT
% =========================================================================

% plot
if plot_comparisons
    
    % downsample for plotting
    ds = 10;
    ms = 12;
    
    % plot values
    figure;
    subplot(1, 2, 1);
    semilogy(x(1:ds:end)*1e3, amp_dec_1D(1:ds:end), 'k.', 'MarkerSize', ms);
    hold on;
    semilogy(x*1e3, amp_dec_1D_an, 'k-');
    semilogy(x(1:ds:end)*1e3, amp_dec_2D(1:ds:end), 'k.', 'MarkerSize', ms);
    semilogy(x*1e3, amp_dec_2D_an, 'k-');
    semilogy(x(1:ds:end)*1e3, amp_dec_3D(1:ds:end), 'k.', 'MarkerSize', ms);
    semilogy(x*1e3, amp_dec_3D_an, 'k-');
    set(gca, 'YLim', [0.01, 10]);
    xlabel('Distance [mm]');
    ylabel('Amplitude [au]');
    title('Decay');
    
    % plot errors
    subplot(1, 2, 2);
    plot(x(ind:end)*1e3, abs(amp_dec_1D(ind:end) - amp_dec_1D_an(ind:end)));
    hold on;
    plot(x(ind:end)*1e3, abs(amp_dec_2D(ind:end) - amp_dec_2D_an(ind:end)));
    plot(x(ind:end)*1e3, abs(amp_dec_3D(ind:end) - amp_dec_3D_an(ind:end)));
    legend('1D', '2D', '3D');
    xlabel('Distance [mm]');
    ylabel('Error [au]');
    title('Error');
    
    scaleFig(1.5, 1);
    
end