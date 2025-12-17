function test_pass = kWaveArray_rotating_sensors(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare the acoustic fields recorded by off-grid
%       sensors with various orientations relative to the grid.
%
% ABOUT:
%       author      - Bradley Treeby and Elliott Wise
%       date        - 16th July 2019
%       last update - 8th June 2021
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby and Elliott Wise

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

% set pass variable
test_pass = true;

% set literals
comparison_thresh = 1;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1;           % [m/s]

% define the source properties
source_freq = 3;                  % source frequency
source_cycles = 3; 
Nrot = 10;                        % number of orientations to simulate
options = {'BLITolerance', 0.01}; % off-grid source options

% computational grid parameters
Nx = 64;
x_size = 4;
Ny = Nx;
Nz = Nx;
dx = x_size/Nx;
dy = dx;
dz = dx;

% =========================================================================
% TWO-DIMENSIONAL SOURCES
% =========================================================================

source_names = {'Line', 'Arc', 'Disc', 'Rect'};

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% get time step
kgrid.makeTime(medium.sound_speed);

% place a point source at the coordinate system's origin
source.p_mask = zeros(Nx, Ny);
source.p_mask(kgrid.x_vec == 0, kgrid.y_vec == 0) = 1;

% create source signal
source.p = toneBurst(1/kgrid.dt, source_freq, source_cycles);

% define orientation angles to test [deg]
angles = linspace(0, 45, Nrot);

% for each source and orientation, compute the error in the sensor data
for source_type = 1:2
    
    % preallocate out matrix
    sensor_data_combined = zeros(Nrot, kgrid.Nt);
    
    % loop through rotations
    for t = 1:Nrot
        
        % define a rotation
        theta = angles(t);
        rot = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
        
        % create empty array
        karray = kWaveArray(options{:});
        
        % compute source grid weights
        switch source_type
            case 1
                
                % set properties for line element
                start_point = rot * [1; -0.5];
                end_point = rot * [1; 0.5];
                
                % add element
                karray.addLineElement(start_point, end_point);
                
            case 2
                
                % set properties for arc element
                arc_pos = rot * [1; 0];
                radius = 0.8;
                diameter = 1;
                focus_pos = rot * [-1; 0];
                
                % add element
                karray.addArcElement(arc_pos, radius, diameter, focus_pos)
                
            case 3
                
                % set properties for disc element
                disc_pos = rot * [1; 0];
                diameter = 1;
                
                % add element
                karray.addDiscElement(disc_pos, diameter)
                
            case 4
                
                % set properties for rect element
                rect_pos = rot * [1; 0];
                Lx = 0.75;
                Ly = 0.75;
           
                % add element
                karray.addRectElement(rect_pos, Lx, Ly, theta)
                
        end

        % assign binary mask from karray to the sensor mask
        sensor.mask = karray.getArrayBinaryMask(kgrid);
        
        % run k-Wave simulation
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
            'PlotSim', plot_simulations);

        % combine data to give one trace per physical array element
        sensor_data_combined(t, :) = karray.combineSensorData(kgrid, sensor_data);       
        
    end
    
    % compute maximum deviations in the sensor data
    deviations = max(100 * (sensor_data_combined - sensor_data_combined(1, :)) ./ max(sensor_data_combined(1, :)), [], 2);
    
    % compute pass
    L_inf = max(deviations(:));
    if (L_inf > comparison_thresh)
        test_pass = false;
    end
    
    % plot the simulated sensor data
    if plot_comparisons
        
        figure;
        plot(sensor_data_combined.');
        xlabel('Time Index');
        ylabel('Pressure');
        title(sprintf('%s source - signals', source_names{source_type}));
        title('Signals');
        
        figure;
        plot((sensor_data_combined - sensor_data_combined(1, :)).');
        xlabel('Time Index');
        ylabel('Pressure');
        title(sprintf('%s source - difference in signals', source_names{source_type}));
        
        figure;
        plot(angles, deviations, '-');
        title(sprintf('%s source - max difference with angle', source_names{source_type}));
        xlabel('Angle [deg]');
        ylabel('[%]');
        
    end
end

% % =========================================================================
% % THREE-DIMENSIONAL SOURCES
% % =========================================================================
% 
% source_names = {'Disc', 'Bowl', 'Rectangle'};
% 
% % create the computational grid
% kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% 
% % record the field at the coordinate system's origin
% sensor_mask = false(Nx, Ny, Nz);
% sensor_mask(kgrid.x_vec == 0, kgrid.y_vec == 0, kgrid.z_vec == 0) = true;
% sensor_data = NaN(1, Nrot);
% 
% % for each source and orientation, compute the error in the sensor data
% for source_type = 1:3
%     
%     % loop through rotations
%     for t = 1:Nrot
%       
%         % define a rotated source base position
%         th = [360 * rand(1), 360 * rand(1), 360 * rand(1)];
%         trans_mat = getAffineMatrix([0, 0, 0], th);
%         base_point = trans_mat(1:3, 1:3) * [1; 0; 0];
%         
%         % compute source grid weights
%         switch source_type
%             case 1
%                 
%                 % set properties for disc element
%                 disc_pos = base_point;
%                 focus_pos = -base_point;
%                 diameter = 1;
%                 
%                 % create array
%                 karray = kWaveArray(options{:});
%                 karray.addDiscElement(disc_pos, diameter, focus_pos)
%                 amp_in = karray.getArrayGridWeights(kgrid);
%                 
%             case 2
%                 
%                 % set properties for bowl element
%                 bowl_pos = base_point;
%                 focus_pos = -base_point;
%                 radius = 1;
%                 diameter = 1;
%                 
%                 % create array
%                 karray = kWaveArray(options{:});
%                 karray.addBowlElement(bowl_pos, radius, diameter, focus_pos)
%                 amp_in = karray.getArrayGridWeights(kgrid);
%                 
%             case 3
%                 
%                 % set properties for rectangular element
%                 rect_pos = base_point;
%                 Lx = 0.75;
%                 Ly = 0.75;
%            
%                 % create array
%                 karray = kWaveArray(options{:});
%                 karray.addRectElement(rect_pos, Lx, Ly, th)
%                 amp_in = karray.getArrayGridWeights(kgrid);
%                 
%         end
%         
%         % set phase to zero
%         phase_in = 0;
%         
%         % compute acoustic field and record sensor_data
%         pressure = acousticFieldPropagator(amp_in, phase_in, dx, freq, sound_speed);
%         sensor_data(t) = pressure(sensor_mask);
%         
%         % plot field
%         if plot_simulations
%             
%             figure(fig);
%             
%             subplot(1, 2, 1);
%             imagesc(abs(pressure(:, :, ceil((Ny + 1)/2))));
%             axis image;
%             colorbar;
%             
%             subplot(1, 2, 2);
%             imagesc(angle(pressure(:, :, ceil((Ny + 1)/2))));
%             axis image;
%             colorbar;
%             
%             colormap(jet(1024));
%             drawnow;
%             
%         end
%         
%     end
%     
%     % compute deviations in the sensor data
%     deviations = 100 * (sensor_data - sensor_data(1)) ./ sensor_data(1);
%     
%     % compute pass
%     L_inf_r = max(real(deviations(:)));
%     L_inf_i = max(imag(deviations(:)));
%     if (L_inf_r > comparison_thresh) || (L_inf_i > comparison_thresh)
%         test_pass = false;
%     end
%     
%     % plot the simulated sensor data
%     if plot_comparisons
%         figure;
%         plot(real(deviations), '-');
%         hold on;
%         plot(imag(deviations), '-');
%         title(sprintf('%s (3D) source', source_names{source_type}));
%         xlabel('Test run');
%         ylabel('Relative error [%]');
%         legend('Real part', 'Imaginary part');
%     end
%     
% end