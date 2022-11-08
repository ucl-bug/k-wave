% Source Smoothing Example
%
% This example illustrates how spatial smoothing can be used to reduce
% discrete sampling artifacts.
%
% For a more detailed discussion of this example, see Treeby, B. E. and
% Cox, B. T., "A k-space Green's function solution for acoustic initial
% value problems in homogeneous media with power law absorption," J.
% Acoust. Soc. Am., vol. 129, no. 6, pp. 3652-3660, 2011.
%
% author: Bradley Treeby
% date: 4th February 2011
% last update: 20th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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
% SOURCE SMOOTHING
% =========================================================================

% create the computational grid
Nx = 256;                   % number of grid points
dx = 0.05e-3;               % grid point spacing [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% define time array
dt = 2e-9;                  % [s]
t_end = 4.26e-6;            % [s]
Nt = round(t_end / dt) + 1;
kgrid.setTime(Nt, dt);

% define the position of the source and sensor
source_sensor_dist = 2e-6;  % [m]
source_pos = Nx/2;
sensor_pos = source_pos + round(medium.sound_speed * source_sensor_dist / dx);

% create a binary sensor mask
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos) = 1;

% set option for symmetric window
symmetric = true;

% create new figure window
figure;

% loop through the window types
for source_index = 1:3
    
    % create delta function source
    source.p0 = zeros(Nx, 1);
    source.p0(source_pos) = 1;
    
    switch source_index
        case 1
            
            % define the window to investigate
            window_type = 'No Window';
            
        case 2

            % define the window to investigate
            window_type = 'Hanning';
            
            % create the filter
            [filt, cg] = getWin(kgrid.Nx, window_type, 'Rotation', true, 'Symmetric', symmetric);
        
            % apply the filter
            source.p0 = real(ifft(fftshift(fftshift(fft(source.p0)) .* filt))) ./ cg;   
            
        case 3

            % define the window to investigate
            window_type = 'Blackman';
            
            % create the filter
            [filt, cg] = getWin(kgrid.Nx, window_type, 'Rotation', true, 'Symmetric', symmetric);
        
            % apply the filter
            source.p0 = real(ifft(fftshift(fftshift(fft(source.p0)) .* filt))) ./ cg;  
            
    end
        
    % run the simulation
    sensor_data = kspaceSecondOrder(kgrid, medium, source, sensor, 'Smooth', false, 'PlotSim', false);
    
    % calculate the amplitude spectrum of the recorded signal
    [f, func_as] = spect(sensor_data, 1/kgrid.dt);
    
    % plot the initial pressure distribution
    subplot(3, 3, (source_index * 3) - 2);
    plot_width = 10;
    stem(0:2*plot_width, source.p0(source_pos - plot_width:source_pos + plot_width), 'k');
    set(gca, 'XLim', [0, 2 * plot_width], 'YLim', [-0.1, 1.1]);
    xlabel('Sample');
    ylabel('Amplitude [au]');
    title([window_type ': Spatial Source Shape']);
    
    % plot the recorded time signals
    subplot(3, 3, (source_index * 3) - 1);
    [~, scale, prefix] = scaleSI(kgrid.t_array(end));
    plot(kgrid.t_array * scale, sensor_data(1, :), 'k-');
    axis tight;
    set(gca, 'YLim', [-0.15, 0.55], 'XLim', [1, 3]);
    ylabel('Amplitude [au]');
    xlabel(['Time [' prefix 's]']);
    title([window_type ': Recorded Time Pulse']);
    
    % plot the amplitude spectrum of the recorded signal
    subplot(3, 3, (source_index * 3));
    [~, scale, prefix] = scaleSI(f(end));
    func_as(1) = 2 * func_as(1);
    plot(f * scale, func_as ./ max(func_as(:)), 'k-');
    set(gca, 'XLim', [0, 20], 'YLim', [0, 1.05]);
    xlabel(['Frequency [' prefix 'Hz]']);
    ylabel('Relative Amplitude Spectrum');
    title([window_type ': Frequency Response']);
    
end

% enlarge figure window
scaleFig(1.75, 2);