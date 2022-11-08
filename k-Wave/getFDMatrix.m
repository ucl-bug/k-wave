function fdm = getFDMatrix(Nx, dx, deriv_order, accuracy_order)
%GETFDMATRIX Create a matrix of finite-difference coefficients.
%
% DESCRIPTION:
%     getFDMatrix returns a matrix of finite-difference coefficients to
%     compute the derivative of the given order with the given numerical
%     accuracy. The coefficients in the center of the grid use centered
%     differences, while the coefficients at the edges of the grid use
%     forward and backward differences. 
%
%     The matrix is generated and returned in sparse format. To convert
%     to a full matrix, use the syntax fdm = full(fdm).
%
%     Example: Compute the gradient of a 1D column vector
%
%         dx = pi/20;
%         x = (0:dx:4*pi).';
%         y = sin(x);
%         fdm = getFDMatrix(length(x), dx);
%         dydx = fdm * y;
%         plot(x, y, 'k-', x, dydx, 'r-');
%
% USAGE:
%     fdm = getFDMatrix(Nx)
%     fdm = getFDMatrix(Nx, dx)
%     fdm = getFDMatrix(Nx, dx, deriv_order, accuracy_order)
%     ...
%
% INPUTS:
%     Nx              - size of data to differentiate
%
% OPTIONAL INPUTS:
%     dx              - spacing between grid-points (default = 1)
%     deriv_order     - order of the derivative to compute, e.g., use 1
%                       to compute df/dx, 2 to compute df^2/dx^2, etc. 
%                       (default = 1)
%     accuracy_order  - order of accuracy for the finite difference
%                       coefficients. Because centered differences are
%                       used, this must be set to an integer multiple of 2
%                       (default = 2)
%
% OUTPUTS:
%     fdm             - matrix of finite difference coefficients
%       
% ABOUT:
%     author          - Bradley Treeby
%     date            - 24th August 2012
%     last update     - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby
%
% See also gradientFD

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

% check if user input for dx is given
if nargin < 2 || isempty(dx)
    dx = 1;
end

% check if user input for deriv_order is given
if nargin < 3 || isempty(deriv_order)
    deriv_order = 1;
elseif ~((deriv_order > 0) && (deriv_order == round(deriv_order)))
    error('Input for deriv_order must be an integer > 0.');
end

% check if user input for accuracy_order is given
if nargin < 4 || isempty(accuracy_order)
    accuracy_order = 2;
end

% check user input for accuracy_order is valid
if rem(accuracy_order, 2)
    error('Input for accuracy_order must be an integer multiple of 2.');
end

% compute the size of the FD stencil width needed for the given accuarcy
% and derivative order
stencil_width = accuracy_order + (deriv_order + rem(deriv_order, 2)) - 1;

% get the finite difference weights
W = getFDWeights(dx, stencil_width, deriv_order);

% calculation the relative indices of the populated diagonals
pop_diag = (0:stencil_width - 1) - (stencil_width - 1)/2;

% generate the FD matrix using central difference coefficients
fdm = spdiags(repmat(W(ceil(end/2), :), Nx, 1), pop_diag, Nx, Nx);

% replace edges with forward and backward difference coefficients
fdm(1:(stencil_width - 1)/2, 1:stencil_width) = W(1:(stencil_width - 1)/2, :);
fdm(end - (stencil_width - 1)/2 + 1:end, end - stencil_width + 1:end) = W(end - (stencil_width - 1)/2 + 1:end, :);


function W = getFDWeights(h, pts, order)
% Compute Finite Difference Weights for a Uniform Grid 
% Function from the MATHWORKS file exchange (ufdwt.m)
%
% Input Parameters:
%
% h     - spacing between FD nodes
% pts   - number of FD points in scheme (3-pt, 5-pt, etc)
% order - order of the derivative operator to compute weights for
%         (note: order<pts-1!)
%         1 computes first derivative differences       
%         2 computes second derivative differences, etc
%  
% Output Parameter:
%
% W is the weight matrix. Each row contains a different set of weights
% (centered or off). If, for example, the number of finite difference points
% is odd, the centered difference weights will appear in the middle row.
%
% Written by: Greg von Winckel - 06/16/04
% Contact: gregvw@chtm.unm.edu

% Modifications:
% 27 August 2012 (Bradley Treeby): moved scaling by h to after computation
% of W to avoid "badly scaled" errors using \ operator

% Copyright (c) 2009, Greg von Winckel
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

N=2*pts-1;  p1=pts-1;

A=repmat((0:p1)',1,N);      B=repmat((-p1:p1),pts,1);

M=(B.^A)./gamma(A+1); 

rhs=zeros(pts,1);   rhs(order+1)=1;

W=zeros(pts,pts);

for k=1:pts
    W(:,k)=M(:,(0:p1)+k)\rhs;
end

W=W';   W(1:pts,:)=W(pts:-1:1,:);

W = W./(h^order);