function sensor_data = kspaceSecondOrder(kgrid, medium, source, sensor, varargin) 
%KSPACESECONDORDER Fast time-domain simulation of wave propagation for homogeneous media.
%
% DESCRIPTION:
%     kspaceSecondOrder simulates the time-domain propagation of linear
%     compressional waves through a one, two, or three dimensional
%     homogeneous acoustic medium given four input structures: kgrid,
%     medium, source, and sensor. The computation is based on an exact
%     second-order k-space model for media with power law absorption. At
%     each time-step (defined by kgrid.dt and kgrid.Nt or kgrid.t_array),
%     the pressure at the positions defined by sensor.mask are recorded and
%     stored. If kgrid.t_array is set to 'auto', this array is
%     automatically generated using the makeTime method of the kWaveGrid
%     class. To prevent wave wrapping, the computational domain can be
%     automatically expanded by a factor of two by setting the optional
%     input 'ExpandGrid' to true.
%
%     An initial pressure distribution can be specified by assigning a
%     matrix (the same size as the computational grid) of arbitrary numeric
%     values to source.p0. An initial pressure gradient can similarly be
%     specified using source.dp0dt. The pressure is returned as an array of
%     time series at the sensor locations defined by sensor.mask. This is
%     specified as a binary matrix (i.e., a matrix of 1's and 0's the same
%     size as the computational grid) representing the grid points within
%     the computational grid that will collect the data. The sensor_data is
%     returned using MATLAB's standard column-wise linear matrix index
%     ordering with the recorded data indexed as
%     sensor_data(sensor_position, time). The final pressure field over the
%     complete computational grid can also be obtained by setting
%     sensor.record to {'p', 'p_final'}. In this case, the output
%     sensor_data is returned as a structure with the outputs appended as
%     the structure fields sensor_data.p and sensor_data.p_final.
%
%     Compared to the first-order simulation functions kspaceFirstOrder1D,
%     kspaceFirstOrder2D, and kspaceFirstOrder3D, kspaceSecondOrder is
%     restricted to homogeneous media and has less functionality. However,
%     it is exact for a homogenous absorbing medium, is more
%     computationally efficient, and allows an initial pressure gradient to
%     be specified.
%
% USAGE:
%     sensor_data = kspaceSecondOrder(kgrid, medium, source, sensor)
%     sensor_data = kspaceSecondOrder(kgrid, medium, source, sensor, ...) 
%
% INPUTS:
%     kgrid               - k-Wave grid object returned by kWaveGrid
%                           containing Cartesian and k-space grid fields 
%     kgrid.t_array       - evenly spaced array of time values [s] (set to
%                           'auto' by kWaveGrid) 
%
%     medium.sound_speed  - homogeneous sound speed within the acoustic
%                           medium [m/s] 
%     medium.alpha_power  - power law absorption exponent
%     medium.alpha_coeff  - power law absorption coefficient 
%                           [dB/(MHz^y cm)]
%     medium.alpha_mode   - optional input to force either the absorption
%                           or dispersion terms in the equation of state to
%                           be excluded; valid inputs are 'no_absorption'
%                           or 'no_dispersion'
%
%     source.p0           - initial pressure within the acoustic medium
%     source.dp0dt        - initial pressure gradient within the acoustic
%                           medium 
%
%     sensor.mask         - binary grid specifying where the pressure is
%                           recorded at each time-step 
%     sensor.record       - cell array of the acoustic parameters to record
%                           in the form sensor.record = {'p'}; valid inputs
%                           are:
%                             'p'       (acoustic pressure)
%                             'p_final' (final pressure field)
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings.
%
%     kgrid               - 
%     'ExpandGrid'        - Boolean controlling whether the grid size is
%                           expanded on two sides to delay the time before
%                           wave wrapping occurs (default = false).  
%     'MeshPlot'          - Boolean controlling whether mesh is used in 
%                           place of imagesc to plot the pressure field
%                           (default = false).  
%     'PlotFrames'        - Boolean controlling whether the pressure field
%                           for each time step is plotted in a new window
%                           (default = false).  
%     'PlotFreq'          - The number of iterations which must pass before
%                           the simulation plot is updated (default = 10). 
%     'PlotScale'         - [min, max] values used to control the scaling
%                           for imagesc (visualisation) (default = [-1, 1]. 
%     'PlotSim'           - Boolean controlling whether the simulation
%                           iterations are progressively plotted (default =
%                           true).  
%     'Smooth'            - Boolean controlling whether source.p0 is
%                           smoothed using smooth before computation
%                           (default = true).
%
% OUTPUTS:
% If sensor.record is not defined by the user:
%     sensor_data         - time varying pressure recorded at the sensor
%                           positions given by sensor.mask 
%
% If sensor.record is defined by the user:
%     sensor_data.p       - time varying pressure recorded at the sensor
%                           positions given by sensor.mask (returned if 'p'
%                           is set)  
%     sensor_data.p_final - final pressure field over the complete domain
%                           (returned if 'p_final' is set) 
%
% ABOUT:
%     author              - Bradley Treeby and Ben Cox
%     date                - 21st August 2008
%     last update         - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2008-2019 Bradley Treeby and Ben Cox
%
% See also kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kWaveGrid, smooth 

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

% start the timer
tic;

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% minimum number of input variables
NUM_REQ_INPUT_VARIABLES = 4;

% optional input defaults
EXPAND_GRID_DEF         = false;
MESH_PLOT_DEF           = false;
PLOT_FRAMES_DEF         = false;
PLOT_FREQ_DEF           = 10;
PLOT_SCALE_DEF          = [-1 1];
PLOT_SIM_DEF            = true;
SMOOTH_P0_DEF           = true;

% create colormap for visualisations
COLOR_MAP               = getColorMap;

% =========================================================================
% CHECK INPUT STRUCTURES
% =========================================================================

% calculate t_array using makeTime if it is not given
if strcmp(kgrid.t_array, 'auto')
    kgrid.makeTime(medium.sound_speed);
end

% check medium fields
checkFieldNames(medium, {'sound_speed', 'density', 'alpha_coeff', 'alpha_power', 'alpha_mode'});
enforceFields(medium, {'sound_speed'});
if isfield(medium, 'density')
    disp('WARNING: medium.density is not a valid input parameter');
end

% check source fields
checkFieldNames(source, {'p0', 'dp0dt'});

% check the absorption mode input is valid
if isfield(medium, 'alpha_mode')
    if ~ischar(medium.alpha_mode) || (~strcmp(medium.alpha_mode, 'no_absorption') && ~strcmp(medium.alpha_mode, 'no_dispersion'))
        error('medium.alpha_mode must be set to ''no_absorption'' or ''no_dispersion''.');
    end
end  

% extract which source fields have been given
source_case = 0;
if isfield(source, 'p0')
    source_case = source_case + 1;
end
if isfield(source, 'dp0dt')
    source_case = source_case + 2;
end

% check sensor fields
checkFieldNames(sensor, {'mask', 'record'});
enforceFields(sensor, {'mask'});

% check the sensor mask is binary
if sum(sensor.mask(:)) ~= numel(sensor.mask) - sum(sensor.mask(:) == 0)
    error('sensor.mask must be a binary grid (numeric values must be 0 or 1).');
end

% set default record options
record.p = true;
record.p_final = false;

% check for sensor.record and set usage flags - if no flags are given, the
% time history of the acoustic pressure is recorded by default
if isfield(sensor, 'record')
    
    % list of allowable record flags
    record.flags = {'p', 'p_final'};

    % check the contents of the cell array are valid inputs
    for record_index = 1:length(sensor.record)
        if ~ismember(sensor.record(record_index), record.flags)
            error(['''' sensor.record{record_index} ''' is not a valid input for sensor.record.']);
        end
    end
    
    % set the usage flags depending on the cell array
    if ismember('p', sensor.record)
        record.p = true;
    else
        % set record.p to false if a user input for sensor.record is given
        % and 'p' is not set 
        record.p = false;
    end
    if ismember('p_final', sensor.record)
        record.p_final = true;
    end 
    
end

% =========================================================================
% CHECK OPTIONAL INPUTS
% =========================================================================

% shorten commonly used field names
t_array = kgrid.t_array;
c = medium.sound_speed;
if isfield(medium, 'alpha_coeff')
    a0 = medium.alpha_coeff;
    y = medium.alpha_power;
else
    a0 = 0;
end 

% assign default input parameters
expand_grid = EXPAND_GRID_DEF;
mesh_plot   = MESH_PLOT_DEF;
plot_frames = PLOT_FRAMES_DEF;
plot_freq   = PLOT_FREQ_DEF;
plot_scale  = PLOT_SCALE_DEF;
plot_sim    = PLOT_SIM_DEF;
smooth_p0   = SMOOTH_P0_DEF;

% replace defaults with user defined values if provided and check inputs    
if nargin < NUM_REQ_INPUT_VARIABLES
    error('Not enough input parameters.');
elseif rem(nargin - NUM_REQ_INPUT_VARIABLES, 2) %**
    error('Optional input parameters must be given as param, value pairs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}                          
            case 'ExpandGrid'
                expand_grid = varargin{input_index + 1};
                if ~islogical(expand_grid)
                    error('Optional input ExpandGrid must be Boolean.');
                end
            case 'MeshPlot'
                mesh_plot = varargin{input_index + 1};
                if ~islogical(mesh_plot)
                    error('Optional input MeshPlot must be Boolean.');
                end
            case 'PlotFrames'
                plot_frames = varargin{input_index + 1};
                if ~islogical(plot_frames)
                    error('Optional input PlotFrames must be Boolean.');
                end                
            case 'PlotFreq'
                plot_freq = varargin{input_index + 1}; 
                if ~(numel(plot_freq) == 1 && isnumeric(plot_freq))
                    error('Optional input PlotFreq must be a single numerical value.');
                end             
            case 'PlotScale'
                plot_scale = varargin{input_index + 1};
                if ~(numel(plot_scale) == 2 && isnumeric(plot_scale))
                    error('Optional input PlotScale must be a 2 element numerical array.');
                end
            case 'PlotSim'
                plot_sim = varargin{input_index + 1};
                if ~islogical(plot_sim)
                    error('Optional input PlotSim must be Boolean.');
                end
            case 'Smooth'
                smooth_p0 = varargin{input_index + 1};
                if ~islogical(smooth_p0)
                    error('Optional input Smooth must be Boolean.');
                end     
            otherwise
                error('Unknown optional input.');
        end
    end
end

% cleanup unused variables
clear *_DEF NUM_REQ_INPUT_VARIABLES time;

% update command line status
disp('Running k-space simulation...'); 

% =========================================================================
% UPDATE COMMAND LINE STATUS
% =========================================================================

% update command line status
disp(['  time steps: ' num2str(length(t_array))]);
switch numDim(kgrid.k)
    case 1
        disp(['  input grid size: ' num2str(kgrid.Nx) ' grid points (' scaleSI(kgrid.x_size) 'm)']);
    case 2
        [~, scale, prefix] = scaleSI(min(kgrid.x_size, kgrid.y_size));
        disp(['  input grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' grid points (' num2str(kgrid.x_size * scale) ' by ' num2str(kgrid.y_size * scale) prefix 'm)']);
    case 3
        [~, scale, prefix] = scaleSI(min([kgrid.x_size, kgrid.y_size, kgrid.z_size]));
        disp(['  input grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' by ' num2str(kgrid.Nz) ' grid points (' num2str(kgrid.x_size * scale) ' by ' num2str(kgrid.y_size * scale) ' by ' num2str(kgrid.z_size * scale) prefix 'm)']); 
end
disp(['  maximum supported frequency: ' scaleSI( kgrid.k_max * min(c(:)) / (2*pi) ) 'Hz']);

% =========================================================================
% CALCULATE FFT OF SOURCE FUNCTIONS AND EXPAND GRID
% =========================================================================

% smooth p0 distribution if required restoring the maximum magnitude
if smooth_p0
    if source_case ~= 2
        disp('  smoothing source.p0 distribution...');  
        source.p0 = smooth(source.p0, true);
    end
    if source_case > 1
        disp('  smoothing source.dp0dt distribution...');  
        source.dp0dt = smooth(source.dp0dt, true);        
    end
end

if ~expand_grid 
    
    % extract wavenumber matrix from kgrid
    k = kgrid.k;
    
    % extract the location of the sensor points
    sensor_index = (sensor.mask == 1);
    
    % compute FFT of source functions
    if source_case ~= 2
        p0_k = fftn(source.p0);
    end
    if source_case > 1
        dp0dt_k = fftn(source.dp0dt);
    end
    
else
    switch numDim(kgrid.k)
        case 1
            
            % create a larger kgrid to prevent wave wrapping
            kgrid_APE = kWaveGrid(kgrid.Nx * 2, kgrid.dx);

            % extract wavenumber matrix from kgrid
            k = kgrid_APE.k;            
            
            % expand the sensor mask
            sensor_mask = [sensor.mask, zeros(size(sensor.mask))];         
            
            % expand the source functions
            if source_case ~= 2
                p0 = zeros(size(k));
                p0(1:kgrid.Nx) = source.p0;
            end
            if source_case > 1
                dp0dt = zeros(size(k));
                dp0dt(1:kgrid.Nx) = source.dp0dt;
            end            
            
        case 2

            % create a larger kgrid to prevent wave wrapping
            kgrid_APE = kWaveGrid(kgrid.Nx * 2, kgrid.dx, kgrid.Ny * 2, kgrid.dy);

            % extract wavenumber matrix from kgrid
            k = kgrid_APE.k;
            
            % expand the sensor mask
            sensor_mask = zeros(size(k));
            sensor_mask(1:kgrid.Nx, 1:kgrid.Ny) = sensor.mask;
            
            % expand the source functions
            if source_case ~= 2
                p0 = zeros(size(k));
                p0(1:kgrid.Nx, 1:kgrid.Ny) = source.p0;
            end
            if source_case > 1
                dp0dt = zeros(size(k));
                dp0dt(1:kgrid.Nx, 1:kgrid.Ny) = source.dp0dt;
            end
            
        case 3
            
            % create a larger kgrid to prevent wave wrapping
            kgrid_APE = kWaveGrid(kgrid.Nx * 2, kgrid.dx, kgrid.Ny * 2, kgrid.dy, kgrid.Nz * 2, kgrid.dz);

            % extract wavenumber matrix from kgrid
            k = kgrid_APE.k;            
            
            % extract the location of the sensor points within the enlarged grid
            sensor_mask = zeros(size(k));
            sensor_mask(1:kgrid.Nx, 1:kgrid.Ny, 1:kgrid.Nz) = sensor.mask;
            
            % expand the source functions
            if source_case ~= 2
                p0 = zeros(size(k));
                p0(1:kgrid.Nx, 1:kgrid.Ny, 1:kgrid.Nz) = source.p0;
            end
            if source_case > 1
                dp0dt = zeros(size(k));
                dp0dt(1:kgrid.Nx, 1:kgrid.Ny, 1:kgrid.Nz) = source.dp0dt;
            end
            
    end
        
    % extract the location of the sensor points within the enlarged grid
    sensor_index = (sensor_mask == 1);     
    
    % compute FFT of source functions
    if source_case ~= 2
        p0_k = fftn(p0);
    end
    if source_case > 1
        dp0dt_k = fftn(dp0dt);
    end
    
    % delete unused variables
    clear kgrid_APE sensor_mask p0 dp0dt source;
    
end

% =========================================================================
% DEFINE TIME PROPAGATORS
% =========================================================================

% shift the wavenumbers
k = ifftshift(k);

% find the index of the zero wave numbers
if source_case > 1
    k_0_index = find(k == 0);
end

% define time propagation function
if a0 ~= 0 
        
    % convert attenuation to nepers
    a0 = db2neper(a0, y);
    
    % define coefficients
    if y == 2 || (isfield(medium, 'alpha_mode') && strcmp('no_dispersion', medium.alpha_mode))
        
        % compute Upsilon
        Upsilon = sqrt(1 - a0.^2 .* c.^(2 .* y) .* k.^(2 .* y - 2));
        
        % reset first value to avoid possible INF when k = 0
        Upsilon(1) = 1;
        
        % precalculate coefficients to speed up computation
        ckU = c .* k .* Upsilon;
        coeff_1 = a0 .* c.^(y + 1) .* k.^y;
        coeff_2 = 0;
        
    elseif (isfield(medium, 'alpha_mode') && strcmp('no_absorption', medium.alpha_mode))
        
        % compute Upsilon
        Upsilon = sqrt(1 - 2 .* a0 .* c.^y .* k.^(y - 1) .* tan(pi .* y ./ 2));
        
        % reset first value to avoid possible INF when k = 0
        Upsilon(1) = 1;
        
        % precalculate coefficients to speed up computation
        ckU = c .* k .* Upsilon;
        coeff_1 = 0;
        coeff_2 = a0 .* c.^y .* k.^(y - 1) ./ Upsilon;
        
        % reset first value to avoid possible INF when k = 0
        coeff_2(1) = 0;        
        
    else
        
        % compute Upsilon
        Upsilon = sqrt(1 - a0.^2 .* c.^(2 .* y) .* k.^(2 .* y - 2) - ...
            2 .* a0 .* c.^y .* k.^(y - 1) .* tan(pi .* y ./ 2));
        
        % reset first value to avoid possible INF when k = 0
        Upsilon(1) = 1;
        
        % precalculate coefficients to speed up computation
        ckU = c .* k .* Upsilon;
        coeff_1 = a0 .* c.^(y + 1) .* k.^y;
        coeff_2 = a0 .* c.^y .* k.^(y - 1) ./ Upsilon;
        
        % reset first value to avoid possible INF when k = 0
        coeff_2(1) = 0;
        
    end
        
else
    
    % precalculate coefficients to speed up computation
    ck = c .* k;
    
end

% =========================================================================
% PREPARE VISUALISATIONS AND STORAGE VARIABLES
% =========================================================================

% preallocate storage variables
if record.p
    sensor_data.p = zeros(sum(sensor.mask(:)), length(t_array));
end

% pre-compute suitable axes scaling factor
if plot_sim
    switch numDim(kgrid.k)
        case 1
            [~, scale, prefix] = scaleSI(max(kgrid.x_vec));
        case 2
            [~, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec])); 
        case 3
            [~, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec; kgrid.z_vec]));
    end
end 

% initialise the figures used for animation if 'PlotSim' is set to 'true'
if plot_sim
    
    % create empty figure
    img = figure;
    
    % create waitbar and shift position so it doesn't overlap the figure window
    pbar = waitbar(0, 'Computing Pressure Field', 'Visible', 'off');
    posn_pbar = get(pbar, 'OuterPosition');
    posn_img = get(img, 'OuterPosition');
    posn_pbar(2) = max(min(posn_pbar(2) - posn_pbar(4), posn_img(2) - posn_pbar(4) - 10), 0);
    set(pbar, 'OuterPosition', posn_pbar, 'Visible', 'on');
    
end

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
disp('  starting time loop...');

% =========================================================================
% LOOP THROUGH TIME STEPS
% =========================================================================

% loop through each value of t
for t_index = 1:length(t_array)

    % extract the time point
    t = t_array(t_index);
    
    % compute pressure field
    switch (source_case + 3 * (a0 ~= 0))
        case 1
            TP_p0 = cos(ck .* t);
            p = ifftn( p0_k .* TP_p0 );  
        case 2
            TP_dp0dt = sin(ck .* t) ./ (ck);
            TP_dp0dt(k_0_index) = t; %#ok<*FNDSB>
            p = ifftn( dp0dt_k .* TP_dp0dt ); 
        case 3
            TP_p0 = cos(ck .* t);
            TP_dp0dt = sin(ck .* t) ./ (ck);
            TP_dp0dt(k_0_index) = t;
            p = ifftn( p0_k .* TP_p0 + dp0dt_k .* TP_dp0dt ); 
        case 4
            TP_p0 = exp(-coeff_1 .* t) .* ( cos(ckU .* t) - coeff_2 .* sin(ckU .* t) );
            p = ifftn( p0_k .* TP_p0 );
        case 5
            TP_dp0dt = sin(ckU .* t) .* exp(-coeff_1 .* t) ./ ckU;
            TP_dp0dt(k_0_index) = t .* exp(-coeff_1(k_0_index)*t);
            p = ifftn( dp0dt_k .* TP_dp0dt );
        case 6
            TP_p0 = exp(-coeff_1*t) .* ( cos(ckU .* t) - coeff_2 .* sin(ckU .* t) );
            TP_dp0dt = sin(ckU .* t) .* exp(-coeff_1 .* t) ./ ckU;
            TP_dp0dt(k_0_index) = t .* exp(-coeff_1(k_0_index) .* t);
            p = ifftn( p0_k .* TP_p0 + dp0dt_k .* TP_dp0dt ); 
    end            

    % extract required data
    if record.p
        sensor_data.p(:, t_index) = real(p(sensor_index));
    end

    % plot data if required
    if plot_frames || ( plot_sim && rem(t_index, plot_freq) == 0 )

        % update progress bar
        waitbar(t_index/length(t_array), pbar);
        drawnow;   

        % update plot
        if plot_frames
            figure;
        end             
        
        p = real(p);

        switch numDim(kgrid.k)
            case 1
                
                % extract the required pressure field from the enlarged grid
                if expand_grid
                    p = p(1:kgrid.Nx);
                end  
                plot(kgrid.x * scale, p);
                xlabel(['x-position [' prefix 'm]']);
                set(gca, 'YLim', plot_scale);
                
            case 2
                
                % extract the required pressure field from the enlarged grid
                if expand_grid
                    p = p(1:kgrid.Nx, 1:kgrid.Ny);
                end                

                if mesh_plot
                    mesh(kgrid.x_vec * scale, kgrid.y_vec * scale, p, 'EdgeColor', 'Black');
                    axis image;
                    set(gca, 'ZLim', plot_scale);
                    axis off;                
                else
                    
                    % add sensor mask onto plot
                    p(sensor.mask == 1) = plot_scale(2);  
                    
                    imagesc(kgrid.y_vec * scale, kgrid.x_vec * scale, p, plot_scale);
                    colormap(COLOR_MAP);
                    axis image;
                    ylabel(['x-position [' prefix 'm]']);
                    xlabel(['y-position [' prefix 'm]']);
                    
                end
          
            case 3
                
                % extract the required pressure field from the enlarged grid
                if expand_grid
                    p = p(1:kgrid.Nx, 1:kgrid.Ny, 1:kgrid.Nz);
                end                  
                planeplot(kgrid, p, '', plot_scale, scale, prefix, COLOR_MAP);
                
        end
    end
end

% =========================================================================
% CLEAN UP
% =========================================================================

% assign output data
if ~isfield(sensor, 'record')
    sensor_data = sensor_data.p;
elseif record.p_final
    switch numDim(kgrid.k)
        case 1
            sensor_data.p_final = p(1:kgrid.Nx);
        case 2
            sensor_data.p_final = p(1:kgrid.Nx, 1:kgrid.Ny);
        case 3
            sensor_data.p_final = p(1:kgrid.Nx, 1:kgrid.Ny, 1:kgrid.Nz);
    end
end

% clean up used figures
if plot_sim
    close(img);
    close(pbar);
end

% update command line status
disp(['  computation completed in ' scaleTime(toc)]);

function planeplot(kgrid, data, data_title, plot_scale, scale, prefix, color_map)
% Subfunction to produce a plot of a three-dimensional matrix through the
% three central planes

subplot(2, 2, 1), imagesc(kgrid.y_vec * scale, kgrid.x_vec * scale, squeeze(data(:, :, kgrid.Nz/2)), plot_scale);
title([data_title 'x-y plane']);
axis image;
subplot(2, 2, 2), imagesc(kgrid.z_vec * scale, kgrid.x_vec * scale, squeeze(data(:, kgrid.Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel(['(All axes in ' prefix 'm)']);
subplot(2, 2, 3), imagesc(kgrid.z_vec * scale, kgrid.y_vec * scale, squeeze(data(kgrid.Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(color_map); 
drawnow;