% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to check
%     and assign optional input parameters (given as 'Param', value pairs)
%     used in the main simulation functions. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 20th November 2014
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2018 Bradley Treeby

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

% replace defaults with user defined values if provided and check inputs    
if num_inputs < NUM_REQ_INPUT_VARIABLES
    error('Not enough input parameters.');
elseif rem(num_inputs - NUM_REQ_INPUT_VARIABLES, 2)
    error('Optional input parameters must be given as param, value pairs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'CartInterp'
                
                % assign input
                cartesian_interp = varargin{input_index + 1};
                
                % check value
                if ~(strcmp(cartesian_interp, 'linear') || strcmp(cartesian_interp, 'nearest'))
                    error('Optional input ''CartInterp'' must be set to ''linear'' or ''nearest''.');
                end
                
            case 'CreateLog'
                
                % assign input
                flags.create_log = varargin{input_index + 1}; 
                
                % check value
                if ~islogical(flags.create_log)
                    error('Optional input ''CreateLog'' must be Boolean.');
                end
                
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
                
                % assign input
                flags.data_recast = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.data_recast)
                    error('Optional input ''DataRecast'' must be Boolean.');
                end
                
            case 'DisplayMask'
                
                % assign input
                display_mask = varargin{input_index + 1};
                if ~(strcmp(display_mask, 'off') || all(size(display_mask) == size(kgrid.k)))
                    error('Optional input ''DisplayMask'' must be the same size as the computational grid or set to ''off''.');
                end
                
                % force mask to be Boolean
                if ~strcmp(display_mask, 'off') 
                    display_mask = (display_mask == 1);
                end
                
            case 'HDFCompressionLevel'
                
                % assign input
                hdf_compression_level = varargin{input_index + 1};
                
                % check value is an integer between 0 and 9
                if (round(hdf_compression_level) ~= hdf_compression_level) ...
                        || (hdf_compression_level > 9) || (hdf_compression_level < 0)
                    error('Optional input ''HDFCompressionLevel'' must be an integer between 0 and 9.');
                end
                    
            case 'LogScale'
                
                % assign input
                flags.plot_scale_log = varargin{input_index + 1};
                
                % check value
                if flags.elastic_code
                    error('Optional input ''LogScale'' is not supported by the elastic code.');
                elseif numel(plot_freq) == 1 && isnumeric(flags.plot_scale_log) && flags.plot_scale_log > 0
                    log_scale_comp_factor = flags.plot_scale_log;
                    flags.plot_scale_log = true;
                elseif ~islogical(flags.plot_scale_log)
                    error('Optional input ''LogScale'' must be Boolean or a single numerical value > 0.');
                end
                
            case 'MeshPlot'
                
                % option only supported in 2D fluid code
                if kgrid.dim == 2 && ~flags.elastic_code
                    flags.mesh_plot = varargin{input_index + 1};
                    if ~islogical(flags.mesh_plot)
                        error('Optional input ''MeshPlot'' must be Boolean.');
                    end
                elseif flags.elastic_code
                    error('Optional input ''MeshPlot'' is not supported by the elastic code.');
                else
                    error('Optional input ''MeshPlot'' only supported in 2D.');
                end
                
            case 'MovieArgs'
                
                % assign input
                movie_args = varargin{input_index + 1};
                
                % check args given as pairs
                if rem(length(movie_args), 2)
                    error('Optional input ''MovieArgs'' must be given as param, value pairs.');
                end
                
            case 'MovieName'
                
                % assign input
                movie_name = varargin{input_index + 1};
                
                % check value
                if ~ischar(movie_name)
                    error('Optional input ''MovieName'' must be a string.');
                end   
                
            case 'MovieType'
                
                % display warning for depracated input
                disp('  WARNING: Optional input ''MovieType'' is no longer supported.');
                
            case 'MovieProfile'
                
                % assign input
                movie_profile = varargin{input_index + 1};
                
            case 'MultiAxialPMLRatio'
                
                % assign input
                multi_axial_PML_ratio = varargin{input_index + 1}; 
                
                % check value
                if ~(numel(multi_axial_PML_ratio) == 1 && isnumeric(multi_axial_PML_ratio) && (multi_axial_PML_ratio >= 0))
                    error('Optional input ''MultiAxialPMLRatio'' must be a single positive value.');
                end
                
            case 'PlotFreq'
                
                % assign input
                plot_freq = varargin{input_index + 1};
                
                % check value
                if ~(numel(plot_freq) == 1 && isnumeric(plot_freq) && (round(plot_freq) == plot_freq) && (plot_freq > 0))
                    error('Optional input ''PlotFreq'' must be a single positive integer value.');
                end
                
            case 'PlotLayout'
                
                % assign input
                flags.plot_layout = varargin{input_index + 1}; 
                
                % check value
                if ~islogical(flags.plot_layout)
                    error('Optional input ''PlotLayout'' must be Boolean.');
                end
                
            case 'PlotPML'
                
                % assign input
                flags.plot_pml = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.plot_pml)
                    error('Optional input ''PlotPML'' must be Boolean.');
                end
                
            case 'PlotScale'
                
                % assign input
                plot_scale = varargin{input_index + 1};
                
                % check value
                if ~flags.elastic_code && ~strcmp(plot_scale, 'auto') && (~(numel(plot_scale) == 2 && isnumeric(plot_scale)))
                    error('Optional input ''PlotScale'' must be a 2 element numerical array or set to ''auto''.');    
                elseif ~strcmp(plot_scale, 'auto') && (~( (numel(plot_scale) == 2 || numel(plot_scale) == 4) && isnumeric(plot_scale)))
                    error('Optional input ''PlotScale'' must be a 2 or 4 element numerical array or set to ''auto''.');    
                end
                
            case 'PlotSim'
                
                % assign input
                flags.plot_sim = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.plot_sim)
                    error('Optional input ''PlotSim'' must be Boolean.');
                end
                
            case 'PMLAlpha'
                
                % check input is correct size
                if length(varargin{input_index + 1}) > kgrid.dim
                    if kgrid.dim > 1
                        error(['Optional input ''PMLAlpha'' must be a 1 or ' kgrid.dim ' element numerical array.']);
                    else
                        error('Optional input ''PMLAlpha'' must be a single numerical value.');
                    end                    
                end
                
                % assign input based on number of dimensions
                switch kgrid.dim
                    case 1
                        pml_x_alpha = varargin{input_index + 1}(1); 
                    case 2
                        pml_x_alpha = varargin{input_index + 1}(1);
                        pml_y_alpha = varargin{input_index + 1}(end);                        
                    case 3
                        pml_x_alpha = varargin{input_index + 1}(1);
                        pml_y_alpha = varargin{input_index + 1}(ceil((end + 1)/2));
                        pml_z_alpha = varargin{input_index + 1}(end);
                end
                
            case 'PMLInside'
                
                % assign input
                flags.pml_inside = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.pml_inside)
                    error('Optional input ''PMLInside'' must be Boolean.');
                end
                
                
            case 'PMLRange'
                
                % assign input
                pml_search_range = varargin{input_index + 1};
                
            case 'PMLSize'
                
                % check for 'auto'
                if ischar(varargin{input_index + 1})
                    if strcmp(varargin{input_index + 1}, 'auto')
                        flags.pml_auto = true;
                    else
                        error(['Optional input ''PMLSize'' must be a 1 or ' kgrid.dim ' element numerical array, or set to ''auto''.']);
                    end
                else
                
                    % check input is correct size
                    if length(varargin{input_index + 1}) > kgrid.dim
                        if kgrid.dim > 1
                            error(['Optional input ''PMLSize'' must be a 1 or ' kgrid.dim ' element numerical array.']);
                        else
                            error('Optional input ''PMLSize'' must be a single numerical value.');
                        end
                    end

                    % assign input based on number of dimensions, rounding to
                    % the nearest integer
                    switch kgrid.dim
                        case 1
                            pml_x_size = round(varargin{input_index + 1}(1));
                        case 2
                            pml_x_size = round(varargin{input_index + 1}(1));
                            pml_y_size = round(varargin{input_index + 1}(end));
                        case 3
                            pml_x_size = round(varargin{input_index + 1}(1));
                            pml_y_size = round(varargin{input_index + 1}(ceil((end + 1)/2)));
                            pml_z_size = round(varargin{input_index + 1}(end));
                    end
                    
                end
                
            case 'RadialSymmetry'
                
                % assign input
                radial_symmetry = varargin{input_index + 1};
                
                % check value
                if ~(strcmp(radial_symmetry, 'WSWA') || strcmp(radial_symmetry, 'WSWS') || strcmp(radial_symmetry, 'WSWA-FFT') || strcmp(radial_symmetry, 'WSWS-FFT'))
                    error('Optional input ''RadialSymmetry'' must be set to ''WSWA'', ''WSWS'', ''WSWA-FFT'', ''WSWS-FFT''.');
                end
                
            case 'RecordMovie'
                
                % assign input
                flags.record_movie = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.record_movie)
                    error('Optional input ''RecordMovie'' must be Boolean.');
                end
                
            case 'StreamToDisk'
                
                % input only supported in 3D fluid code
                if ~flags.elastic_code && kgrid.dim == 3
                    
                    % assign input
                    flags.stream_to_disk = varargin{input_index + 1};
                    
                    % check value
                    if ~(numel(flags.stream_to_disk) == 1 && ( isnumeric(flags.stream_to_disk) || islogical(flags.stream_to_disk) ))
                        error('Optional input ''StreamToDisk'' must be a single scalar or Boolean value.');
                    end
                    
                    % if given as a Boolean, replace with the default
                    % number of time steps
                    if islogical(flags.stream_to_disk) && (flags.stream_to_disk ~= false)
                        flags.stream_to_disk = STREAM_TO_DISK_STEPS_DEF;
                    end
                    
                else
                    error('Optional input ''StreamToDisk'' is currently only compatible with 3D fluid simulations.');
                end
                
            case 'SaveToDisk'
                
                % not supported in 1D
                if kgrid.dim ~= 1
                    
                    % check for h5create function
                    if ~exist('h5create', 'file')
                        error('Optional input ''SaveToDisk'' requires the function h5create which is included with MATLAB 2011a or later.');
                    end
                    
                    % assign and check input
                    flags.save_to_disk = varargin{input_index + 1};
                    if ~(islogical(flags.save_to_disk) || ischar(flags.save_to_disk))
                        error('Optional input ''SaveToDisk'' must be Boolean or a String.');
                    end
                    
                else
                    error('Optional input ''SaveToDisk'' is not compatible with 1D simulations.');
                end
                
                % assign input
                if islogical(flags.save_to_disk) && flags.save_to_disk
                    flags.save_to_disk = SAVE_TO_DISK_FILENAME_DEF;
                end
                
            case 'SaveToDiskExit'
                
                % not supported in 1D
                if kgrid.dim ~= 1
                    
                    % assign input
                    flags.save_to_disk_exit = varargin{input_index + 1};
                    
                    % check value
                    if ~islogical(flags.save_to_disk_exit)
                        error('Optional input ''SaveToDiskExit'' must be Boolean.');
                    end
                    
                else
                    error('Optional input ''SaveToDisk'' is not compatible with 1D simulations.');
                end
                
            case 'ScaleSourceTerms'
                
                % assign input
                flags.scale_source_terms = varargin{input_index + 1};
                
                % check value
                if ~islogical(flags.scale_source_terms)
                    error('Optional input ''ScaleSourceTerms'' must be Boolean.');
                end
                
            case 'Smooth'
                
                % check input size
                if length(varargin{input_index + 1}) > 3 || ~islogical(varargin{input_index + 1})
                    error('Optional input ''Smooth'' must be a 1, 2 or 3 element Boolean array.');
                end
                
                % assign input
                flags.smooth_p0 = varargin{input_index + 1}(1);
                flags.smooth_c0 = varargin{input_index + 1}(ceil((end + 1)/2));
                flags.smooth_rho0 = varargin{input_index + 1}(end);
                
            case 'UseFD'
                
                % input only supported in 1D fluid code
                if kgrid.dim == 1 && ~flags.elastic_code
                    
                    % assign input
                    flags.use_finite_difference = varargin{input_index + 1};
                    
                    % check value
                    if ~((islogical(flags.use_finite_difference) && ~flags.use_finite_difference) || ...
                            (isnumeric(flags.use_finite_difference) && (flags.use_finite_difference == 2 || flags.use_finite_difference == 4)))
                        error('Optional input ''UseFD'' must be set to 2, 4, or false.');
                    end
                    
                else
                    error('Optional input ''UseFD'' only supported in 1D.');
                end  
                
            case 'UsekSpace'
                
                % assign input
                flags.use_kspace = varargin{input_index + 1}; 
                
                % check value
                if ~islogical(flags.use_kspace)
                    error('Optional input ''UsekSpace'' must be Boolean.');
                end
                
            case 'UseSG'
                
                % assign input
                flags.use_sg = varargin{input_index + 1}; 
                
                % check value
                if ~islogical(flags.use_sg)
                    error('Optional input ''UseSG'' must be Boolean.');
                end
                
            otherwise
                error(['Unknown optional input ' varargin{input_index} '.']);
        end
    end
end

% automatically assign the PML size to give small prime factors
if flags.pml_auto
    if flags.pml_inside
        error('''PMLSize'' set to ''auto'' is only supported with ''PMLInside'' set to false.');
    else
        
        % get optimal pml size
        if flags.axisymmetric
            pml_size_temp = getOptimalPMLSize(kgrid, pml_search_range, radial_symmetry(1:4));
        else
            pml_size_temp = getOptimalPMLSize(kgrid, pml_search_range);
        end
        
        % assign to individual variables
        switch kgrid.dim
            case 1
                pml_x_size = pml_size_temp(1);
            case 2
                pml_x_size = pml_size_temp(1);
                pml_y_size = pml_size_temp(2);
            case 3
                pml_x_size = pml_size_temp(1);
                pml_y_size = pml_size_temp(2);
                pml_z_size = pml_size_temp(3);
        end
        
        % cleanup unused variables
        clear pml_size_temp;        
        
    end
end