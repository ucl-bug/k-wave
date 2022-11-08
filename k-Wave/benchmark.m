function benchmark(datacast)
%BENCHMARK Run performance benchmark.
%
% DESCRIPTION:
%     benchmark performs an analysis of the time taken to run simulations
%     using kspaceFirstOrder3D with increasing grid sizes. The default
%     simulation uses a heterogeneous absorbing medium, a binary sensor
%     mask with 100 points, 1000 time steps, 'PlotSim' set to false, and
%     grid sizes varying from 32^3 to 256^3. The times are computed from 3
%     averages and the function stops when memory errors are encountered or
%     the benchmarking is complete. The computational times and settings
%     are stored in a MATLAB formatted binary file using save. 
%
%     Note, before running, the MATLAB workspace is cleared and any open
%     MATLAB windows are closed. 
%
% INPUTS:
%     datacast    - option for 'DataCast' (default = 'off')
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 27th October 2010
%     last update - 25th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby
%
% See also kspaceFirstOrder3D

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

%#ok<*ASGLU>
%#ok<*UNRCH>

% define literals
USE_DIARY       = false;
DATA_CAST_DEF   = 'off';

% close windows
close all hidden;

% store the matlab and k-Wave versions
options = getComputerInfo;
options.report_mem_usage = strncmp(computer, 'PCWIN', 5);

% set the data cast option, e.g.,
%   'off'             double precision
%   'single'          single precision
%   'gpuArray-double' double precision on GPU using Parallel Computing Toolbox
%   'gpuArray-single' single precision on GPU using Parallel Computing Toolbox
if nargin == 0
    options.data_cast = DATA_CAST_DEF;
else
    options.data_cast = datacast;
end

% simulation options
options.heterogeneous_media     = true;
options.absorbing_media         = true;
options.nonlinear_media         = false;
options.binary_sensor_mask      = true;
options.smooth                  = [false, false, false];
options.number_sensor_points    = 100;
options.plot_sim                = false;
options.number_time_points      = 1000;
options.num_averages            = 3;

% set the remaining input options
input_args = {'PlotSim', options.plot_sim, 'PMLSize', 10, 'DataCast', options.data_cast, 'Smooth', options.smooth}; 

% set the computational grid dimensions to use
start_size = 32;
x_scale_array = [1 2 2 2 4 4 4 8 8 8 16 16];
y_scale_array = [1 1 2 2 2 4 4 4 8 8 8  16];
z_scale_array = [1 1 1 2 2 2 4 4 4 8 8  8];

% create filenames for the diary and benchmark data
filename = ['benchmark_data-' options.computer_name '-' options.data_cast '-' options.date];
error_reached = false;

% start diary if required
if USE_DIARY
    diaryname = ['benchmark_diary-' options.computer_name '-' options.data_cast '-' options.date]; 
    diary(diaryname);   
end

% ambient medium properties
c0      = 1500;   
rho0    = 1000;
a0      = 0.75;
y       = 1.5;
BonA    = 6;

% run simulations until completion or a memory error is thrown
try
    for size_index = 1:length(x_scale_array)
        
        % extract the scale parameters
        xscale = x_scale_array(size_index);
        yscale = y_scale_array(size_index);
        zscale = z_scale_array(size_index);

        % scale parameter for source radius
        scale = min([xscale yscale zscale]);

        % create the computational grid
        Nx = start_size * xscale;    
        Ny = start_size * yscale;    
        Nz = start_size * zscale;    
        dx = 22e-3 / Nx;
        dy = 22e-3 / Ny;
        dz = 22e-3 / Nz;       
        kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

        % define the properties of the propagation medium
        if options.heterogeneous_media
            medium.sound_speed = c0*ones(Nx, Ny, Nz);     
            medium.sound_speed(1:Nx/4, :, :) = c0 * 1.2;

            medium.density = rho0*ones(Nx, Ny, Nz); 
            medium.density(:, Ny/4:end, :) = rho0 * 1.2;
        else
            medium.sound_speed = c0;
            medium.density = rho0;
        end
        if options.absorbing_media
            medium.alpha_coeff = a0;
            medium.alpha_power = y;
        end
        if options.nonlinear_media
            medium.BonA = BonA;
        end

        % create initial pressure distribution using makeBall
        ball_magnitude = 10;
        ball_x_pos = Nx/2;  
        ball_y_pos = Ny/2;  
        ball_z_pos = Nz/2;  
        ball_radius = 2 * scale;  
        source.p0 = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);
        
        % smooth the initial pressure distribution if not done within the
        % simulation function
        if ~options.smooth(1)
            source.p0 = smooth(source.p0, true);
        end
        
        % define sensor points to collect the data
        sensor.mask = makeCartSphere(10e-3, options.number_sensor_points);
        
        if options.binary_sensor_mask
            sensor.mask = cart2grid(kgrid, sensor.mask);
        end
        
        % define the time array
        kgrid.makeTime(max(medium.sound_speed(:))); 
        kgrid.setTime(options.number_time_points, kgrid.dt);
        
        % reset the timer count
        loop_time = 0;
        
        % reset memory count
        loop_mem_usage = 0;
        
        for loop_num = 1:options.num_averages
            
            % run the simulation 
            if options.report_mem_usage
                
                % start the timer
                tic;

                % run simulation and collect memory usage
                [sensor_data, mem] = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
                
                % stop the timer
                elapsed_time = toc;
                
                % compute the average memory used by MATLAB during the simulation
                loop_mem_usage = (loop_mem_usage * (loop_num - 1) + mem.user.MemUsedMATLAB)/loop_num;
            
                % store the memory usage
                mem_usage(size_index) = loop_mem_usage;  %#ok<AGROW,NASGU>
                
            else
                
                % start the timer
                tic;
                
                % run simulation
                kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});  
                
                % stop the timer
                elapsed_time = toc;
                
            end
            
            % compute the average loop time
            loop_time = (loop_time * (loop_num - 1) + elapsed_time) / loop_num;

            % store the loop time
            comp_time(size_index) = loop_time; %#ok<AGROW,NASGU>

            % store the computation size
            comp_size(size_index) = Nx * Ny * Nz;  %#ok<AGROW,NASGU>

            % save the loop variables each loop in case of memory errors
            if options.report_mem_usage
                save(filename, 'comp_size', 'comp_time', 'mem_usage', 'options');
            else
                save(filename, 'comp_size', 'comp_time', 'options');
            end
            
        end
    end
catch ME
    disp(' ');
    disp('Memory limit reached or error encountered, exiting benchmark. Error message:');
    error_reached = true;
end

% switch off diary
if USE_DIARY
    diary off;
end

% display the error
if error_reached
    disp([ '  ' ME.message]);
    disp(' ');
end