function varargout = angularSpectrum(input_plane, dx, dt, z_pos, medium, varargin)
% ANGULARSPECTRUM Project time-domain input plane using the angular spectrum method.
%
% DESCRIPTION:
%     angularSpectrum projects a 2D input plane (given as a 3D matrix of
%     time series at each spatial position) using the angular spectrum
%     method. The time series are decomposed into spectral components and
%     then each frequency is propagated using the spectral propagator with
%     angular restriction described in reference [1]. 
%
%     The time signals are computed in retarded time, where the time series
%     for each plane are offset by z/c0 and the input plane corresponds to
%     z = 0. The projections are calculated directly from the input plane
%     (no spatial stepping is used).
%
%     For linear projections in a lossless medium, just the sound speed can
%     be specified. For projections in a lossy medium, the parameters are
%     given as fields to the input structure medium. 
%
%     Two datasets can be returned. The first, pressure_max, contains a 3D
%     matrix of the maximum (temporal peak) pressure across the 2D planes
%     specified by z_pos, and is indexed as (x_ind, y_ind, plane_index).
%     The second, pressure_time, contains the time-varying pressure signals 
%     across the 2D planes specified by z_pos, and is indexed as (x_ind,
%     y_ind, t_ind, plane_index).
%
%     To compute the maximum pressure field over an isotropic domain with
%     Nz grid points (assuming the source plane is aligned with z_ind = 1)
%     without storing the time-series (which can be memory consuming), use
%     the syntax:
%
%           pressure_max = angularSpectrum(input_plane, dx, dt, (0:(Nz - 1)) * dx, c0)
%
%     To compute the time series over another parallel plane separated by a
%     distance of z_pos, use the syntax
%
%           [~, pressure_time] = angularSpectrum(input_plane, dx, dt, z_pos, c0)
%
% USAGE:
%     pressure_max = angularSpectrum(input_plane, dx, dt, z_pos, c0)
%     pressure_max = angularSpectrum(input_plane, dx, dt, z_pos, medium, ...)
%     [pressure_max, pressure_time] = angularSpectrum(input_plane, dx, dt, z_pos, c0)
%     [pressure_max, pressure_time] = angularSpectrum(input_plane, dx, dt, z_pos, medium, ...)
%     ...
%
% INPUTS:
%     input_plane          - 3D matrix containing the time varying pressure
%                            over a 2D input plane indexed as (x, y, t)
%                            [Pa].
%     dx                   - Spatial step between grid points in the input
%                            plane [m].
%     dt                   - Temporal step between time points in the input
%                            plane [s].
%     z_pos                - Vector specifying the relative z-position of
%                            the planes to which the data is projected [m].
%
%     c0                   - Medium sound speed [m/s].
%             OR
%     medium.sound_speed   - Medium sound speed [m/s].
%     medium.alpha_power   - Power law absorption exponent.
%     medium.alpha_coeff   - Power law absorption coefficient
%                            [dB/(MHz^y cm)].
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'AngularRestriction' - Boolean controlling whether angular
%                            restriction is used as described in [1]
%                            (default = true).
%     'DataCast'           - String input of the data type that variables
%                            are cast to before computation. For example,
%                            setting to 'single' will speed up the
%                            computation time (due to the improved
%                            efficiency of fft2 and ifft2 for this data
%                            type). This variable is also useful for
%                            utilising GPU parallelisation the Parallel
%                            Computing Toolbox by setting 'DataCast' to
%                            'gpuArray-single' (default = 'off').
%     'DataRecast'         - Boolean controlling whether the output data
%                            is cast back to double precision. If set to
%                            false, sensor_data will be returned in the
%                            data format set using the 'DataCast' option
%                            (default = false).
%     'FFTLength'          - Length of the FFT used to compute the angular
%                            spectrum (default = 1 + the next power of two
%                            larger than the grid size).
%     'GridExpansion'      - Grid padding used to increase the accuracy of
%                            the projection. The grid expansion is removed
%                            before returning the calculated pressure to
%                            the user (default = 0).
%     'Plot'               - Boolean controlling whether the field is
%                            plotted at each z step (default = true if
%                            z_pos contains more than one value). 
%     'Reverse'            - Boolean controlling whether the projection is
%                            in the forward (false) or backward (true)
%                            direction (default = false).
%
% OUTPUTS:
%     pressure_max         - 3D matrix of maximum pressure (temporal peak) 
%                            across the 2D planes specified by z_pos,
%                            indexed as (x_ind, y_ind, plane_index) [Pa].
%     pressure_time        - 4D matrix of the time varying pressure signals
%                            across the 2D planes specified by z_pos,
%                            indexed as (x_ind, y_ind, t_ind, plane_index)
%                            [Pa].
% 
% ABOUT:
%     author               - Bradley Treeby
%     date                 - 27th February 2018
%     last update          - 19th February 2019
%
% REFERENCES:
%     [1] Zeng, X., & McGough, R. J. (2008). Evaluation of the angular
%     spectrum approach for simulations of near-field pressures. The
%     Journal of the Acoustical Society of America, 123(1), 68-76.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby
%
% See also angularSpectrumCW

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

% start timer
start_time = clock;

% =========================================================================
% INPUT CHECKING
% =========================================================================

% define defaults
angular_restriction = true;
grid_expansion      = 0;
fft_length          = 'auto';
data_cast           = 'off';
data_recast         = false;
reverse_proj        = false;
absorbing           = false;
plot_updates        = true;
loops_for_time_est  = 5;
record_time_series  = false;

% check for the number of outputs
if nargout == 2
    record_time_series = true;
end

% replace with user defined values if provided
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'AngularRestriction'
                angular_restriction = logical(varargin{input_index + 1});
            case 'DataCast'
                
                % assign input
                data_cast = varargin{input_index + 1};
                
                % check list of valid inputs
                if ~ischar(data_cast)
                    error('Optional input ''DataCast'' must be a string.');
                elseif ~(strcmp(data_cast, 'off') || strcmp(data_cast, 'double') ...
                        || strcmp(data_cast, 'single') || strcmp(data_cast, 'gpuArray-single') ... 
                        || strcmp(data_cast, 'gpuArray-double'))
                    error('Invalid input for ''DataCast''.');
                end
                
                % replace double with off
                if strcmp(data_cast, 'double')
                    data_cast = 'off';
                end
                
                % create empty string to hold extra cast variable for use
                % with the parallel computing toolbox
                data_cast_prepend = '';
                
                % replace PCT options with gpuArray
                if strcmp(data_cast, 'gpuArray-single')
                    data_cast = 'gpuArray';
                    data_cast_prepend = 'single';
                elseif strcmp(data_cast, 'gpuArray-double')
                    data_cast = 'gpuArray';
                end
                
                if strcmp(data_cast, 'gpuArray')
                    
                    % check the PCT is installed and the version is 2012a or
                    % later (verLessThan only works in release 7.4 or later)
                    v = ver;
                    if verLessThan('matlab', '7.14') || ~ismember('Parallel Computing Toolbox', {v.Name})
                        error('The Parallel Computing Toolbox for MATLAB 2012a or later is required for ''DataCast'' set to ''gpuArray-single'' or ''gpuArray-double''.');
                    end

                    % cleanup unused variables
                    clear v;
                    
                end
                
            case 'DataRecast'
                data_recast = logical(varargin{input_index + 1});
            case 'FFTLength'
                fft_length = round(varargin{input_index + 1});
            case 'GridExpansion'
                grid_expansion = round(varargin{input_index + 1});
            case 'Plot'
                plot_updates = logical(varargin{input_index + 1});
            case 'Reverse'
                reverse_proj = logical(varargin{input_index + 1});
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check for structured medium input
if isstruct(medium)
    
    % force the sound speed to be defined
    if ~isfield(medium, 'sound_speed')
        error('medium.sound_speed must be defined when specifying medium properties using a structure.');
    end
    
    % assign the sound speed
    c0 = medium.sound_speed;
    
    % assign the absorption
    if isfield(medium, 'alpha_coeff') || isfield(medium, 'alpha_power')
        
        % enforce both absorption parameters
        if ~(isfield(medium, 'alpha_coeff') && isfield(medium, 'alpha_power'))
            error('Both medium.alpha_coeff and medium.alpha_power must be defined for an absorbing medium.');
        end
        
        % assign flag
        absorbing = true;
        
    end
    
else
    
    % assign the sound speed
    c0 = medium;
    
end

% check time step is sufficient
if (c0 * dt / dx) > 1
    error('Maximum supported frequency in temporal sampling is lower than maximum supported frequency in spatial sample (CFL > 1).');
end

% =========================================================================
% PRE-PROCESSING
% =========================================================================

% get grid size
[Nx, Ny, Nt] = size(input_plane);
Nz = length(z_pos);

% turn off plotting if only one value of z_pos
if numel(z_pos) == 1
    plot_updates = false;
end

% get scale factor for grid size
[~, scale, prefix] = scaleSI(min(Nx * dx, Ny * dx));

% update command line status
disp('Running angular spectrum projection...');
disp(['  start time: ' datestr(start_time)]);
disp(['  input plane size: ' num2str(Nx) ' by ' num2str(Ny) ' grid points (' num2str(scale * Nx * dx) ' by ' num2str(scale * Ny * dx) prefix 'm)']);
disp(['  grid expansion: ' num2str(grid_expansion) ' grid points']);

% reverse time signals if stepping backwards
if reverse_proj
    input_plane = flip(input_plane, 3);
end

% expand input
if grid_expansion > 0
    input_plane = expandMatrix(input_plane, [grid_expansion, grid_expansion, 0], 0);
    [Nx, Ny, Nt] = size(input_plane);
end

% get FFT size
if ischar(fft_length) && strcmp(fft_length, 'auto')
    fft_length = 2.^(nextpow2(max([Nx, Ny])) + 1);
end

% update command line status
disp(['  FFT size: ' num2str(fft_length) ' points']);
disp(['  maximum supported frequency: ' scaleSI( c0 / (2 * dx) ) 'Hz']);
disp(['  input signal length: ' num2str(Nt) ' time points (' scaleSI(Nt * dt) 's)']);

% create wavenumber vector
N = fft_length;
if rem(N, 2) == 0
    k_vec = ((-N/2):(N/2-1)) .* 2 * pi ./ (N * dx);
else
    k_vec = (-(N-1)/2:(N-1)/2) .* 2 * pi ./ (N * dx);
end

% force middle value to be zero in case 1/Nx is a recurring
% number and the series doesn't give exactly zero
k_vec(floor(N/2) + 1) = 0;

% shift wavenumbers to be in the correct order for FFTW
k_vec = ifftshift(k_vec);

% create wavenumber grids
[kx, ky] = meshgrid(k_vec, k_vec);
 
% precompute term for angular restriction
sqrt_kx2_ky2 = sqrt(kx.^2 + ky.^2);

% preallocate maximum pressure output
pressure_max = zeros(Nx, Ny, Nz);

% preallocate time series output
if record_time_series
    pressure_time = zeros(Nx, Ny, Nt, Nz);
end

% =========================================================================
% DECOMPOSE INPUT INTO FREQUENCY COMPONENTS
% =========================================================================

% compute the fft
input_plane_w_fft = fft(input_plane, [], 3);

% reduce to a single sided spectrum where the number of unique points for
% even numbered FFT lengths is given by N/2 + 1, and for odd (N + 1)/2
num_unique_pts = ceil((Nt + 1) / 2);
input_plane_w_fft = input_plane_w_fft(:, :, 1:num_unique_pts);

% create the frequency axis variable
f_vec = (0:size(input_plane_w_fft, 3) - 1) / (dt * Nt);

% compute frequencies to propagate
f_vec_prop = f_vec(f_vec < (c0 / (2 * dx)));

% preallocate loop variable
pressure_time_step = zeros(Nx, Ny, length(f_vec));

% =========================================================================
% DATA CASTING
% =========================================================================

% cast variables to the output type
if ~strcmp(data_cast, 'off')
    
    % update command line status
    disp(['  casting variables to ' data_cast ' type...']);
    
    % list of variables to cast
    cast_variables = {'kx', 'ky', 'z_pos', 'input_plane_w_fft', 'pressure_max', 'pressure_time_step'};
    
    % additional variable if storing snapshots
    if record_time_series
        cast_variables = [cast_variables, {'pressure_time'}];
    end
    
    % additional variables used if absorbing
    if absorbing
        cast_variables = [cast_variables, {'k'}];
    end
        
    % additional variables used for angular restriction
    if angular_restriction
        cast_variables = [cast_variables, {'sqrt_kx2_ky2', 'fft_length', 'dx'}];
    end
    
    % loop through, and change data type
    for cast_index = 1:length(cast_variables)
        eval([cast_variables{cast_index} ' = ' data_cast '(' data_cast_prepend '(' cast_variables{cast_index} '));']);
    end
    
end

% =========================================================================
% Z-LOOP
% =========================================================================

% open figure window
if plot_updates
    sim_fig = figure;
    scaleFig(1.5, 1);
end

% update command line status
loop_start_time = clock;
disp(['  precomputation completed in ' scaleTime(etime(loop_start_time, start_time))]);
disp('  starting z-step loop...');

% loop over z-positions
for z_index = 1:Nz
    
    % get current z value
    z = z_pos(z_index);
    
    % if set to zero, just store the input plane
    if z == 0
        
        % store maximum pressure
        pressure_max(:, :, z_index) = max(input_plane, [], 3);
        
        % store time series data if required
        if record_time_series
            pressure_time(:, :, :, z_index) = input_plane;
        end
        
    else

        % loop over frequencies
        for f_index = 1:length(f_vec_prop)

            % compute wavenumber at driving frequency
            k = 2 * pi * f_vec(f_index) / c0;

            % compute wavenumber grid
            kz = sqrt(k.^2 - complex(kx.^2 + ky.^2));

            % compute spectral propagator (Eq. 6)
            H = conj(exp(1i .* z .* kz));

            % account for attenuation (Eq. 11)
            if absorbing
                
                % convert attenuation to Np/m
                alpha_Np = db2neper(medium.alpha_coeff, medium.alpha_power) * (2 * pi * f_vec(f_index))^medium.alpha_power;
                
                % apply attenuation to propagator
                if alpha_Np ~= 0
                    H = H .* exp(-alpha_Np .* z .* k ./ kz);
                end
                
            end

            % apply angular restriction
            if angular_restriction

                % size of computational domain [m]
                D = (fft_length - 1) * dx;

                % compute angular threshold (Eq. 10)
                kc = k * sqrt(0.5 * D.^2 ./ (0.5 * D.^2 + z.^2));

                % apply threshold to propagator
                H(sqrt_kx2_ky2 > kc) = 0;

            end

            % compute forward Fourier transform of input plane
            input_plane_xy_fft = fft2(input_plane_w_fft(:, :, f_index), fft_length, fft_length);

            % compute phase shift for retarded time
            ret_time = exp(1i * 2 * pi * f_vec(f_index) * z / c0);

            % compute projected field
            pressure_step = ifft2(input_plane_xy_fft .* H, fft_length, fft_length);
            pressure_time_step(:, :, f_index) = pressure_step(1:Nx, 1:Ny) .* ret_time;

        end

        % form double sided amplitude spectrum in correct order for FFTW
        % (FFT of real data is conjugate symmetric) 
        if rem(Nt, 2)
            pressure_time_step_exp = cat(3, pressure_time_step, flip(conj(pressure_time_step(:, :, 2:end)), 3));
        else
            pressure_time_step_exp = cat(3, pressure_time_step, flip(conj(pressure_time_step(:, :, 2:end - 1)), 3));
        end

        % take inverse Fourier tranform to recover time domain data
        pressure_time_step_exp = real(ifft(pressure_time_step_exp, [], 3));

        % store maximum pressure
        pressure_max(:, :, z_index) = max(pressure_time_step_exp, [], 3);

        % store time series data if required
        if record_time_series
            pressure_time(:, :, :, z_index) = pressure_time_step_exp;
        end

        % update command line status
        if z_index == loops_for_time_est
            disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time) * Nz / z_index) '...']);
        end
        
        % plot updates
        if plot_updates
            figure(sim_fig);
            subplot(1, 2, 1);
            imagesc(reshape(pressure_time_step_exp, [], Nt));
            axis square;
            xlabel('time index');
            ylabel('linear grid index');
            colorbar;
            title('Time Series');
            subplot(1, 2, 2);
            imagesc(max(pressure_time_step_exp, [], 3));
            axis image;
            xlabel('x-grid index');
            ylabel('y-grid index');
            title(['Maximum Pressure (z-step ' num2str(z_index) ' of ' num2str(Nz) ')']);
            colorbar;
            drawnow;
        end        
        
    end
    
end

% update command line status
disp(['  simulation completed in ' scaleTime(etime(clock, loop_start_time))]);

% =========================================================================
% POST PROCESSING
% =========================================================================

% close plot
if plot_updates
    close(sim_fig);
end

% trim grid expansion
if grid_expansion > 0
    pressure_max = pressure_max(1 + grid_expansion:end - grid_expansion, 1 + grid_expansion:end - grid_expansion, :);
    if record_time_series
        pressure_time = pressure_time(1 + grid_expansion:end - grid_expansion, 1 + grid_expansion:end - grid_expansion, :, :);
    end
end

% reverse time signals and grid if stepping backwards
if reverse_proj
    pressure_max = flip(pressure_max, 3);
    if record_time_series
        pressure_time = flip(pressure_time, 3);
    end
end

% cast output back to double precision
if data_recast
    
    % update command line status
    disp('  recasting variables to double...');
    
    % recast data
    if strncmp(data_cast, 'gpuArray', 8)
        pressure_max = double(gather(pressure_max));
        if record_time_series
            pressure_time = double(gather(pressure_time));
        end
    else
        pressure_max = double(pressure_max);
        if record_time_series
            pressure_time = double(pressure_time);
        end
    end
end

% assign outputs
varargout{1} = pressure_max;
if record_time_series
    varargout{2} = pressure_time;
end

% update command line status
disp(['  total computation time ' scaleTime(etime(clock, start_time))]);