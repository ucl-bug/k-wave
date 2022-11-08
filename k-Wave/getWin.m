function [win, cg] = getWin(N, type, varargin)
%GETWIN Return a frequency domain windowing function.
%
% DESCRIPTION:
%     getWin returns a 1D, 2D, or 3D frequency domain window of the
%     specified type of the given dimensions. By default, higher 
%     dimensional windows are created using the outer product. The windows
%     can alternatively be created using rotation by setting the optional
%     input 'Rotation' to true. The coherent gain of the window can also be
%     returned. 
%
% USAGE:
%     win = getWin(N, type)
%     win = getWin(N, type, ...)
%     [win, cg] = getWin(N, type)
%     [win, cg] = getWin(N, type, ...)
%
% INPUTS:
%     N           - number of samples, use
%
%                   N = Nx for 1D
%                   N = [Nx Ny] for 2D
%                   N = [Nx Ny Nz] for 3D
%
%     type        - window type. Supported values are
%
%                   'Bartlett'
%                   'Bartlett-Hanning'   
%                   'Blackman'              
%                   'Blackman-Harris'
%                   'Blackman-Nuttall'
%                   'Cosine'
%                   'Flattop'
%                   'Gaussian'
%                   'HalfBand'
%                   'Hamming'     
%                   'Hanning'  
%                   'Kaiser' 
%                   'Lanczos'
%                   'Nuttall'            
%                   'Rectangular'                
%                   'Triangular'
%                   'Tukey'
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Plot'      - Boolean controlling whether the window is displayed
%                   (default = false).
%
%     'Param' -     Control parameter for the Tukey, Blackman, Gaussian,
%                   and Kaiser windows:
%
%                   Tukey: taper ratio (default = 0.5)
%                   Blackman: alpha (default = 0.16)
%                   Gaussian: standard deviation (default = 0.5)
%                   Kaiser: alpha (default = 3)
%
%     'Rotation'  - Boolean controlling whether 2D and 3D windows are
%                   created via rotation or the outer product (default =
%                   false). Windows created via rotation will have edge
%                   values outside the window radius set to the first
%                   window value. 
%
%     'Symmetric' - Boolean controlling whether the window is symmetrical
%                   (default = true). If set to false, a window of length N
%                   + 1 is created and the first N points are returned. For
%                   2D and 3D windows, 'Symmetric' can be defined as a
%                   vector defining the symmetry in each matrix dimension.
%
%     'Square'    - Boolean controlling whether the window is forced to
%                   be square (default = false). If set to true and Nx
%                   and Nz are not equal, the window is created using the
%                   smaller variable, and then padded with zeros. 
%
% OUTPUTS:
%     win         - the window
%     cg          - the coherent gain of the window
%
% ABOUT:
%     author      - Bradley E. Treeby
%     date        - 10th March 2010
%     last update - 25th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby
%
% See also spectrum

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

% set usage defaults
num_req_input_variables = 2;
plot_win = false;
rotation = false;
symmetric = true;
square = false;
ylim = [0, 1];
if strcmp(type, 'Tukey')
    param = 0.5;
    param_ub = 1; 
    param_lb = 0;
elseif strcmp(type, 'Blackman')
    param = 0.16;
    param_ub = 1; 
    param_lb = 0;    
elseif strcmp(type, 'Gaussian')
    param = 0.5;
    param_ub = 0.5; 
    param_lb = 0;   
elseif strcmp(type, 'Kaiser')
    param = 3;
    param_ub = 100; 
    param_lb = 0;      
else
    param = 0;
end

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Plot'
                plot_win = varargin{input_index + 1};
            case 'Param'
                param = varargin{input_index + 1}(1);
                if param > param_ub
                    param = param_ub;
                elseif param < param_lb
                    param = param_lb;
                end
            case 'Rotation'
                rotation = varargin{input_index + 1};
                if ~islogical(rotation)
                    error('Optional input Rotation must be Boolean.');
                end                  
            case 'Symmetric'
                
                % assign input
                symmetric = varargin{input_index + 1};
                
                % check type
                if ~islogical(symmetric)
                    error('Optional input Symmetric must be Boolean.');
                end   
                
                % check size
                if ~((numel(symmetric) == 1) || (numel(symmetric) == numel(N)))
                    error('Optional input Symmetric must have 1 or numel(N) elements.');
                end
                
            case 'Square'
                square = varargin{input_index + 1};
                if ~islogical(square)
                    error('Optional input Square must be Boolean.');
                end               
            otherwise
                error('Unknown optional input.');
        end
    end
end

% set any required input options for recursive function calls
if ismember(type, {'Tukey', 'Blackman', 'Kaiser', 'Gaussian'})
    input_options = {'Param', param};
else
    input_options = {};
end

% if a non-symmetrical window is required, enlarge the window size (note,
% this expands each dimension individually if symmetric is a vector)
N = N + 1 .* ~symmetric;

% if a square window is required, replace grid sizes with smallest size and
% store a copy of the original size
if square && (length(N) ~= 1)
    N_orig = N;
    L = min(N);
    N(:) = L;
end

% create the window
switch length(N)
    case 1
        n = 0:N-1;
        switch type
            case 'Bartlett'
                win = (2 / (N - 1) * ( (N - 1) / 2 - abs(n - (N - 1) / 2))).';
            case 'Bartlett-Hanning'
                win = (0.62 - 0.48 * abs(n / (N - 1) - 1/2) - 0.38 * cos(2 * pi * n / (N - 1))).';     
            case 'Blackman'              
                win = cosineSeries(n, N, [(1 - param)/2, 0.5, param/2]);
            case 'Blackman-Harris'
                win = cosineSeries(n, N, [0.35875, 0.48829, 0.14128, 0.01168]);
            case 'Blackman-Nuttall'
                win = cosineSeries(n, N, [0.3635819, 0.4891775, 0.1365995, 0.0106411]);
            case 'Cosine'
                win = (cos(pi * n / (N - 1) - pi/2)).';   
            case 'Flattop'
                win = cosineSeries(n, N, [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368]);
                ylim = [-0.2, 1];
            case 'Gaussian'
                win = (exp(-0.5 * ( (n - (N - 1) / 2) / (param * (N - 1) / 2) ).^2)).'; 
            case 'HalfBand'
                win = ones(N, 1);
                ramp_length = round(N / 4);
                ramp = 1/2 + 9/16 * cos(pi * (1:ramp_length) ./ (2 * ramp_length)) - 1/16 *cos(3 * pi * (1:ramp_length) ./ (2 * ramp_length));               
                win(1:ramp_length) = fliplr(ramp);
                win(end - ramp_length + 1:end) = ramp;
            case 'Hamming'
                win = (0.54 - 0.46 * cos(2 * pi * n / (N - 1))).';
            case 'Hanning'
                win = (0.5 - 0.5 * cos(2 * pi * n / (N - 1))).';
            case 'Kaiser'
                win = (besseli(0, pi * param * sqrt(1 - (2 * n / (N - 1) - 1).^2 )) / besseli(0, pi * param)).';  
            case 'Lanczos'
                win = sinc(2 * pi * n / (N - 1) - pi).';
            case 'Nuttall'
                win = cosineSeries(n, N, [0.3635819, 0.4891775, 0.1365995, 0.0106411]);             
            case 'Rectangular'
                win = ones(N, 1);
            case 'Triangular'
                win = (2 / N * (N / 2 - abs(n - (N - 1) / 2))).';
            case 'Tukey'
                win = ones(N, 1);
                index = 0:(N - 1) * param / 2;
                param = param * N;
                win(1:length(index)) = 0.5 * (1 + cos(2 * pi / param * (index - param / 2)));              
                win(end:-1:end - length(index) + 1) = win(1:length(index));                
            otherwise
                error(['Unknown window type: ' type]);
        end
        
        % trim the window if required
        N = N - 1.* ~symmetric;
        win = win(1:N);
        
        % calculate the coherent gain
        cg = sum(win)/N;  
        
    case 2

        % create the 2D window
        if rotation
            
            % create the window in one dimension using getWin recursively
            L = max(N);
            win_lin = getWin(L, type, input_options{:});
            
            % create the reference axis
            radius = (L - 1)/2;
            ll = linspace(-radius, radius, L);
            
            % create the 2D window using rotation
            xx = linspace(-radius, radius, N(1));
            yy = linspace(-radius, radius, N(2));
            [x, y] = ndgrid(xx, yy);
            r = sqrt(x.^2 + y.^2);
            r(r > radius) = radius;
            win = interp1(ll, win_lin, r); 
            win(r <= radius) = interp1(ll, win_lin, r(r <= radius));
                        
        else
            
            % create the window in each dimension using getWin recursively
            win_x = getWin(N(1), type, input_options{:});            
            win_y = getWin(N(2), type, input_options{:});
            
            % create the 2D window using the outer product
            win = (win_y * win_x').';
           
        end
        
        % trim the window if required
        N = N - 1.* ~symmetric;
        win = win(1:N(1), 1:N(2));
        
        % calculate the coherent gain
        cg = sum(win(:))/prod(N);        
        
    case 3
        
        % create the 3D window
        if rotation
            
            % create the window in one dimension using getWin recursively
            L = max(N);
            win_lin = getWin(L, type, input_options{:});
            
            % create the reference axis
            radius = (L - 1)/2;
            ll = linspace(-radius, radius, L);
            
            % create the 3D window using rotation
            xx = linspace(-radius, radius, N(1));
            yy = linspace(-radius, radius, N(2));
            zz = linspace(-radius, radius, N(3));
            [x, y, z] = ndgrid(xx, yy, zz);
            r = sqrt(x.^2 + y.^2 + z.^2);
            r(r > radius) = radius;
            win = interp1(ll, win_lin, r); 
            win(r <= radius) = interp1(ll, win_lin, r(r <= radius));            
            
        else        
        
            % create the window in each dimension using getWin recursively
            win_x = getWin(N(1), type, input_options{:});
            win_y = getWin(N(2), type, input_options{:});
            win_z = getWin(N(3), type, input_options{:});

            % create the 2D window using the outer product
            win_2D = (win_x * win_z');

            % create the 3D window
            win = zeros(N(1), N(2), N(3));
            for index = 1:N(2)
                win(:, index, :) = win_2D(:, :) * win_y(index);
            end
            
        end
        
        % trim the window if required
        N = N - 1.* ~symmetric;
        win = win(1:N(1), 1:N(2), 1:N(3));
           
        % calculate the coherent gain
        cg = sum(win(:)) / prod(N);          
        
    otherwise
        error('Invalid input for N, only 1-, 2-, and 3-D windows are supported.');
end

% enlarge the window if required
if square && (length(N) ~= 1)
    L = N(1);
    win_sq = win;
    win = zeros(N_orig);
    switch length(N)
        case 2
            index1 = round((N(1) - L)/2) + 1;
            index2 = round((N(2) - L)/2) + 1;
            win(index1:index1 + L - 1, index2:index2 + L - 1) = win_sq;
        case 3
            index1 = floor((N_orig(1) - L)/2) + 1;
            index2 = floor((N_orig(2) - L)/2) + 1;
            index3 = floor((N_orig(3) - L)/2) + 1;
            win(index1:index1 + L - 1, index2:index2 + L - 1, index3:index3 + L - 1) = win_sq;
    end
end

if plot_win
    switch length(N)
        case 1
            figure;
            plot(win, 'k-');
            xlabel('Samples');
            ylabel('Amplitude');
            box on;
            grid on;
            axis tight;
            set(gca, 'YLim', ylim);
            title([num2str(N) ' sample ' type ' window']);
        case 2
            figure;
            surf(win);
            xlabel('y samples');
            ylabel('x samples');
            zlabel('Amplitude');
            axis tight;
            set(gca, 'DataAspectRatio', [1 1 1/max(N)], 'ZLim', ylim);
            title([num2str(N(1)) ' \times ' num2str(N(2)) ' sample ' type ' window']);
        case 3
            figure;
            
            subplot(2, 2, 1);
            imagesc(squeeze(win(:, :, round(N(3)/2))));
            axis image;
            xlabel('y');
            ylabel('x');
            title('x-y slice');
            
            subplot(2, 2, 2);
            imagesc(squeeze(win(:, round(N(2)/2), :)));
            axis image;
            title('x-z slice');
            xlabel('z');
            ylabel('x');
            
            subplot(2, 2, 3)
            imagesc(squeeze(win(round(N(1)/2), :, :)));
            axis image;
            title('y-z slice');
            xlabel('z');
            ylabel('y');
            
            if exist('sgtitle') %#ok<EXIST>
            	sgtitle([num2str(N(1)) ' \times ' num2str(N(2)) ' \times ' num2str(N(3)) ' sample ' type ' window']);
            end
    end
end

function series = cosineSeries(n, N, coeffs)
% Sub-function to calculate a summed filter cosine series.

series = coeffs(1);
for index = 2:length(coeffs)
    series = series + (-1)^(index-1) * coeffs(index) * cos((index - 1) * 2 * pi * n / (N - 1));
end
series = series.';