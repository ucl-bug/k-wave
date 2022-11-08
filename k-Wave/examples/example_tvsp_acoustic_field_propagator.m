% Simulating CW Fields Using The Acoustic Field Propagator Example
%
% This example demonstrates how to use acousticFieldPropagator to simulate
% the steady-state pressure field from a steered line array in 2D without
% time stepping.
%
% For a more detailed discussion of this example and the underlying
% techniques, see B. E. Treeby, J. Budisky, E. S. Wise, J. Jaros, and B. T.
% Cox, "Rapid calculation of acoustic fields from arbitrary continuous-wave
% sources," The Journal of the Acoustical Society of America, vol. 143, no.
% 1, pp.529-537, 2018.
%
% author: Bradley Treeby
% date: 4th April 2019
% last update: 4th April 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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
% SIMULATION
% =========================================================================

% define the computational grid
Nx = 192;           % number of grid points in the x (row) direction
Ny = 192;           % number of grid points in the y (column) direction
dx = 0.1e-3;    	% grid point spacing [m]

% define source frequency and sound speed
f0 = 1e6;           % source frequency [Hz]
c0 = 1500;          % sound speed [m/s]

% allocate empty amplitute and phase matrices
amp_in = zeros(Nx, Ny);
phase_in = zeros(Nx, Ny);

% define source aperture as a line array at the top of the grid
aperture_width = 60;
y1 = Ny/2 - aperture_width/2 + 1;
y2 = Ny/2 + aperture_width/2;
x1 = 1;

% assign constant amplitude across the line array
amp_in(x1, y1:y2) = 1;

% define angle sweep from -60 to 60 degrees and back
angle_array = [-60:2:60, 58:-2:-60];

% calculate position of each grid point in the array
el_pos = dx * (0:aperture_width - 1);

% open new figure window
figure;
scaleFig(1.5, 1);

% loop through angles
for index = 1:length(angle_array)

    % get the current steering angle
    steering_angle = angle_array(index);
    
    % calculate phase offset for each grid point in the line array based on
    % element position and steering angle, and assign to the line array
    phase_in(x1, y1:y2) = 2 * pi * f0 * el_pos * sind(steering_angle) / c0;

    % compute beam pattern
    [amp_out, phase_out] = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0);
    
    % plot amplitude
    subplot(1, 2, 1);
    imagesc(amp_out);
    axis image;
    title(['Amplitude (angle = ' num2str(angle_array(index)) ' deg)']);
    colorbar
    axis off;

    % plot phase
    subplot(1, 2, 2);
    imagesc(phase_out);
    axis image;
    title('Phase');
    colorbar
    axis off;
    
    colormap(parula(256));
    drawnow;
    
end