function I_out = vesselFilter(I_in, grid_spacing, scales, varargin)
%VESSELFILTER Frangi's 3D vessel filter.
% 
% DESCRIPTION:
%     vesselFilter filters a 3D image using Frangi's vessel filtering
%     algorithm [1-3]. The algorithm works by calculating the Hessian
%     matrix (containing second order gradients) at each image voxel. The
%     eigenvalues of this matrix are then ordered and used to classify
%     whether the voxel is part of a vessel. The algorithm is performed
%     over multiple scales by convolving the input image with a Gaussian of
%     the given input scales. The final output is taken as the maximum of
%     the vessel filtered image across all scales.
%
% USAGE:
%     I_out = vesselFilter(I_in, grid_spacing, scales)
%     I_out = vesselFilter(I_in, grid_spacing, scales, ...)
%
% INPUTS:
%     I_in         - 3D image data to filter
%     grid_spacing - [dx, dy, dz] in [m]
%     scales       - array of scales to use [m]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'alpha'      - Value of sensitivity parameter for metric that
%                    distinguishes between plate-like and other
%                    structures (vessel-like or ball-like) (default = 0.5).
%     'beta'       - Value of sensitivity parameter for metric that
%                    distinguishes between ball-like and other structures
%                    (vessel-like or plate-like) (default = 0.5).
%     'c'          - Value used to scale the sensitivity parameter for
%                    the noise metric (the sensitivity parameter itself
%                    is calculated automatically based on the magnitudes
%                    of the eigenvalues) (default = 1).
%     'gamma '     - normalisation factor for scale-space derivatives
%                    (default = 1)
%     'DisplayUpdates'
%                  - Boolean controlling whether command line updates are
%                    displayed (default = true).
%     'Plot'       - Boolean controlling whether a maximum intensity
%                    projection of the vessel filtered image at each
%                    scale is displayed (default = false)
%     'Colormap'   - colormap to use if 'Plot' is set to true (default =
%                    flipud(gray))
%       
% OUTPUTS:
%     I_out        - vessel filtered image
%
% ABOUT:
%     author       - Bradley Treeby and Tanmayi Oruganti
%     date         - 9th May 2012
%     last update  - 7th June 2017
%
%     This function is based on original code by Dean Barratt and Yipeng
%     Hu, CMIC, UCL, 2006-2009.
%
% REFERENCES:
%     [1] A. F. Frangi, W. J. Niessen, K. L. Vincken, M. A. Viergever 
%     (1998) "Multiscale vessel enhancement filtering," MICCAI, pp. 130-7.
%     [2] R. Manniesing, M.A. Viergever, W.J. Niessen (2006) "Vessel
%     enhancing diffusion: A scale space representation of vessel
%     structures," Med. Image Anal. 10, pp. 815-25. 
%     [3] T. Oruganti, J. Laufer, B. E. Treeby (2013) "Vessel filtering of
%     photoacoustic images," Proc. of SPIE, vol. 8581, p. 85811W-1.
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby and Tanmayi Oruganti
%
% See also interpftn, smooth

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

% set default parameters
num_req_input_variables = 3;    % minimum number of input variables
alpha = 0.5;                    % value of sensitivity parameter for plate vs vessel/ball measure
beta = 0.5;                     % value of sensitivity parameter for ball vs vessel/plate measure
c_scale = 1;                    % scale value of sensitivity parameter for noise
gamma = 0;                      % normalisation factor for scale-space derivatives
origin_smoothness_const = 1e-6; % constant for Manniesing's origin smoothness term
disp_updates = true;            % Boolean controlling display of command line updates
plot_updates = false;           % Boolean controlling MIP display 
color_map = flipud(gray);       % colormap to use if images are displayed 
cardan_tol = 0.001;             % tolerance for finding complex roots using cardanRoots

% check user defined inputs
% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'alpha'
                alpha = varargin{input_index + 1};
            case 'beta'             
                beta = varargin{input_index + 1};
            case 'c'
                c_scale = varargin{input_index + 1};
            case 'gamma'
                gamma = varargin{input_index + 1};
            case 'Display'
                disp_updates = varargin{input_index + 1};
            case 'DisplayUpdates'
                disp_updates = varargin{input_index + 1};                
            case 'Plot'
                plot_updates = varargin{input_index + 1};
            case 'Colormap'
                color_map = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% get size of input volume
sz = size(I_in);

% create k-space grid to use for spectral derivatives
kgrid = kWaveGrid(sz(1), grid_spacing(1), sz(2), grid_spacing(2), sz(3), grid_spacing(3));

% force wavenumber vectors to be in the correct directions for bsxfun
kx_vec = kgrid.kx_vec;
ky_vec = kgrid.ky_vec.';
kz_vec = permute(kgrid.kz_vec, [2 3 1]);

% compute 3D fft of input image
Ik = fftn(I_in);

% preallocate empty output image
I_out = zeros(size(I_in));

% loop through the scales defined by the user
for scale_index = 1:length(scales)
    
    % extract the current scale
    scale = scales(scale_index);
    
    % update command line status
    if disp_updates
        disp(['Computing Vesselness filter for scale ' num2str(scale) ':']);
        tic, fprintf('  Calculating derivatives... ');
    end

    % create normalised frequency domain gaussian at current scale
    FG = exp( (-0.5*scale^2)*(kgrid.kx.^2 + kgrid.ky.^2 + kgrid.kz.^2) );

    % compute components of Hessian matrix convolved with Gaussian
    Hxx = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * kx_vec).^2, FG)) )), [], 1);
    Hyy = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * ky_vec).^2, FG)) )), [], 1);
    Hzz = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * kz_vec).^2, FG)) )), [], 1);
    Hxy = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * kx_vec), bsxfun(@times, (1i * ky_vec), FG))) )), [], 1);
    Hxz = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * kx_vec), bsxfun(@times, (1i * kz_vec), FG))) )), [], 1);
    Hyz = reshape(real(ifftn( Ik .* ifftshift( bsxfun(@times, (1i * ky_vec), bsxfun(@times, (1i * kz_vec), FG))) )), [], 1);

    % update command line status
    if disp_updates
        toc, tic, fprintf('  Computing eigenvalues... ');
    end

    % compute trace of Hessian matrix
    P2 = -(Hxx + Hyy + Hzz);
    
    % compute principal minors of Hessian matrix
    M11 = Hyy .* Hzz - Hyz.^2;
    M22 = Hzz .* Hxx - Hxz.^2;
    M33 = Hxx .* Hyy - Hxy.^2;
    P1 = (M11 + M22 + M33);
    
    % compute determinant of Hessian matrix
    P0 = - Hxx .* M11 ...
         + Hxy .* (Hxy .* Hzz - Hyz .* Hxz) ...
         - Hxz .* (Hxy .* Hyz - Hyy .* Hxz);

    % find eigenvalues using roots of 3rd order polynomial
    P3 = 1;
    eigenval_mat = cardanRoots(P3, P2, P1, P0, cardan_tol).';    
    
    % update command line status
    if disp_updates
        toc, tic, fprintf('  Sorting eigenvalues... ');
    end

    % sort eigenvalues based on their absolute values
    [~, ind] = sort(abs(eigenval_mat));
    ind = bsxfun(@plus, (0:numel(I_in) - 1) * 3, ind);
    eigenval_mat = eigenval_mat(ind);

    % update command line status
    if disp_updates
        toc, tic, fprintf('  Evaluating vesselness measure... ');
    end

    % calculate noise parameter dynamically if not given
    Rs = sqrt( abs(eigenval_mat(1, :)).^2 + abs(eigenval_mat(2, :)).^2 + abs(eigenval_mat(3, :)).^2 );
    c = c_scale * 0.5 * max(Rs(:)); 
    
    % calculate overall vesselness measure and scale by Lindeberg?s constant
    %   term 1: plate vs vessel/ball measure
    %   term 2: ball vs vessel/plate measure
    %   term 3: noise measure (Frobenius matrix norm of the Hessian)
    %   term 4: Manniesing's origin smoothness term
    V = scale.^(gamma) .* ...
        (1 - exp(-(abs(eigenval_mat(2, :)) ./ abs(eigenval_mat(3, :))).^2 ./ (2 * alpha^2))) .* ...
        (exp(-(abs(eigenval_mat(1, :)) ./ sqrt(abs(eigenval_mat(2, :) .* eigenval_mat(3, :)))).^2 ./ (2 * beta^2))) .* ...
        (1 - exp(-Rs.^2 ./ (2 * c^2))) .* ...
        (exp(-(2 * origin_smoothness_const.^2) ./ (abs(eigenval_mat(2, :)) .* eigenval_mat(3, :).^2)));
    
    % eliminate structures with positive eigenvalues as these indicate
    % regions of low optical absorption on a background of high optical
    % absorption 
    V(eigenval_mat(2, :) > 0 | eigenval_mat(3, :) > 0) = 0;

    % reshape to be the correct size
    V = reshape(V, sz);
    
    % take maximum response over multiple scales
    I_out = max(I_out, V);
    
    % display time elapsed
    if disp_updates
        toc
    end
    
    % plot steps if required
    if plot_updates
        
        % get suitable axis scaling
        [~, axis_scale, axis_prefix] = scaleSI(max([sz(1) * grid_spacing(1), sz(2) * grid_spacing(2)]));
        
        % plot figure
        figure;
        imagesc(kgrid.y_vec * axis_scale, kgrid.x_vec * axis_scale,  max(V, [], 3));
        axis image;
        colormap(color_map);
        title(['Vessel Filtered Image (MIP), Scale = ' num2str(scale)]);
        ylabel(['x [' axis_prefix 'm]']);
        xlabel(['y [' axis_prefix 'm]']);
        colorbar;
        drawnow;
        
        % plot final multi-scale image
        if scale_index == length(scales)            
            figure;
            imagesc(kgrid.y_vec * axis_scale, kgrid.x_vec * axis_scale,  max(I_out, [], 3));
            axis image;
            colormap(color_map);
            title('Vessel Filtered Image (MIP), Multi-Scale');
            ylabel(['x [' axis_prefix 'm]']);
            xlabel(['y [' axis_prefix 'm]']);
            colorbar;
            drawnow;
        end
    end
end

% display total time elapsed
if disp_updates
    disp(['  Total computation time ' scaleTime(etime(clock, start_time))]);
end

% end of function
end

% =========================================================================
% SUB FUNCTIONS
% =========================================================================

function roots = cardanRoots(varargin)
% Find roots of third order polynomials using Cardan's formula
%
% INPUT
%   P: (n x 4) array, each row corresponds to coefficients of each
%   polynomial, P(:,1)*x^3 + P(:,2)*x^2 + P(:,3)*x + P(:,4)
% OUTPUT
%   roots: (n x 3) array, each row correspond to the roots of P
%
% To adjust the parameter below which the the discriminant is considerered
% as nil, use
%   CardanRoots(P, tol)
% Adjusting tol is useful to avoid the real roots become complex due to
% numerical accuracy. The default TOL is 0
%
% http://www.sosmath.com/algebra/factor/fac11/fac11.html
% http://mathforum.org/dr.math/faq/faq.cubic.equations.html
%
% See also: roots, ParabolaRoots, eig3
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 20-May-2010

% Copyright (c) 2010, Bruno Luong
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Adjustable parameter
tol = 0;

if nargin<4
    P = varargin{1};
    a = P(:,1);
    b = P(:,2);
    c = P(:,3);
    d = P(:,4);
    if nargin>=2
        tol = varargin{2};
    end
else
    [a, b, c, d] = deal(varargin{1:4});
    if nargin>=5
        tol = varargin{5};
    end
end

if ~isequal(a,1)
    b = b./a;
    c = c./a;
    d = d./a;
end

b2 = b.^2;

p = -b2/3 + c;
q = ((2/27)*b2-(1/3)*c).*b + d;
delta = q.^2 + (4/27)*p.^3;

% Three cases of discriminant sign
iscmplx = imag(p) | imag(q);
notcmplx = ~iscmplx;
deltanull = notcmplx & abs(delta)<tol; % = 0
deltaneg = notcmplx & delta<0;
deltapos = notcmplx & ~(deltanull | deltaneg);

n = size(delta,1);
roots = zeros(n, 3, class(delta));

if any(deltanull)
    idx = find(deltanull);
    roots(idx,:) = CardanNull(p(idx), q(idx));
end

if any(deltaneg)
    idx = find(deltaneg);
    roots(idx,:) = CardanNeg(q(idx), delta(idx));
end

if any(deltapos)
    idx = find(deltapos);
    roots(idx,:) = CardanPos(q(idx), delta(idx));
end

if any(iscmplx)
    idx = find(iscmplx);
    roots(idx,:) = CardanCmplx(p(idx), q(idx), delta(idx));
end

roots = bsxfun(@minus, roots, b/3);

end


function roots = CardanNull(p, q)

    S1 = 3*q./p;
    S2 = -0.5*S1;
    roots = [S1 S2 S2];  % double real solutions

end


function roots = CardanNeg(q, delta)

    alfa = -q;
    beta = sqrt(-delta);
    r2 = alfa.^2-delta;
    rho = (4^(1/3))*exp(log(r2)/6);
    theta = atan2(beta,alfa)/3;
    alfa = rho.*cos(theta);
    beta = rho.*sin(theta);
    S1 = alfa;
    x = (-0.5)*alfa;
    y = (sqrt(3)/2)*beta;
    S2 = x-y;
    S3 = x+y;
    roots = [S1 S2 S3];

end


function roots = CardanPos(q, delta)

    sqrtdelta = sqrt(delta);
    u3 = (-q+sqrtdelta)/2;
    v3 = (-q-sqrtdelta)/2;

    % Cubic roots of u3 and v3
    u = sign(u3).*exp(log(abs(u3))/3);
    v = sign(v3).*exp(log(abs(v3))/3);

    S1 = u+v;
    % Complex solutions
    j = complex(-0.5,sqrt(3)/2);
    j2 = complex(-0.5,-sqrt(3)/2);
    S2 = j*u+j2*v;
    S3 = conj(S2);
    roots = [S1 S2 S3];

end


function roots = CardanCmplx(p, q, delta)

    sqrtdelta = sqrt(delta);
    u3 = (-q+sqrtdelta)/2;
    v3 = (-q-sqrtdelta)/2;

    % we need u*v = -p/3
    p = (-1/3)*p;
    iu = abs(u3)>abs(v3);
    u = zeros(size(u3),class(u3));
    v = zeros(size(v3),class(v3));
    u(iu) = exp(log(u3(iu))/3);
    v(iu) = p(iu)./u(iu);
    v(~iu) = exp(log(v3(~iu))/3);
    u(~iu) = p(~iu)./v(~iu);

    S1 = u+v;

    j = complex(-0.5,sqrt(3)/2);
    j2 = complex(-0.5,-sqrt(3)/2);
    S2 = j*u+j2*v;
    S3 = j2*u+j*v;
    roots = [S1 S2 S3];

end