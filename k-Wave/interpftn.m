function y = interpftn(x, sz, win)
%INTERPFTN Resample data using Fourier interpolation.
%
% DESCRIPTION:
%     interpftn resamples an N-D matrix to the size given in sz using
%     Fourier interpolation.
%
% USAGE:
%     y = interpftn(x, sz)
%     y = interpftn(x, sz, win)
%
% INPUTS:
%     x           - matrix to interpolate
%     sz          - new size
%
% OPTIONAL INPUTS:
%     win         - name of windowing function to use
%
% OUTPUTS:
%     y           - resampled matrix       
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 20th October 2010
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby
%
% See also interpft

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

% check inputs for window
if nargin == 2
    win = 0;
end

% extract the size of the input matrix
x_sz = size(x);

% check enough coefficients have been given
if sum(x_sz ~= 1) ~= numel(sz)
    error('The number of scaling coefficients must equal the number of dimensions in x.');
end

% compute the interpolation over each matrix dimension using the inbuilt
% interpft (dimensions with no interpolation required are skipped)
y = x;
for p = 1:numel(sz)
    if sz(p) ~= x_sz(p)
        y = interpft(y, sz(p), p, win);
    end
end

function y = interpft(x, ny, dim, win)
%INTERPFT 1-D interpolation using FFT method.
%   Y = INTERPFT(X,N) returns a vector Y of length N obtained
%   by interpolation in the Fourier transform of X. 
%
%   If X is a matrix, interpolation is done on each column.
%   If X is an array, interpolation is performed along the first
%   non-singleton dimension.
%
%   INTERPFT(X,N,DIM) performs the interpolation along the
%   dimension DIM.
%
%   Assume x(t) is a periodic function of t with period p, sampled
%   at equally spaced points, X(i) = x(T(i)) where T(i) = (i-1)*p/M,
%   i = 1:M, M = length(X).  Then y(t) is another periodic function
%   with the same period and Y(j) = y(T(j)) where T(j) = (j-1)*p/N,
%   j = 1:N, N = length(Y).  If N is an integer multiple of M,
%   then Y(1:N/M:N) = X.
%
%   Example: 
%      % Set up a triangle-like signal signal to be interpolated 
%      y  = [0:.5:2 1.5:-.5:-2 -1.5:.5:0]; % equally spaced
%      factor = 5; % Interpolate by a factor of 5
%      m  = length(y)*factor;
%      x  = 1:factor:m;
%      xi = 1:m;
%      yi = interpft(y,m);
%      plot(x,y,'o',xi,yi,'*')
%      legend('Original data','Interpolated data')
%
%   Class support for data input x:
%      float: double, single
%  
%   See also INTERP1.

%   Robert Piche, Tampere University of Technology, 10/93.
%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.15.4.4 $  $Date: 2006/12/15 19:27:56 $

% y = zeros(ny, 1);
% y(1:length(x)) = fftshift(fft(x));
% y = real(ifft(ifftshift(y)));

% Bug fix, additional documentation, and win option
% author        - Bradley Treeby
% date          - 20th October 2010
% last update   - 1st April 2011

% if nargin==2,
%     
%     % this forces x to be a column vector, nshifts is so it can be changed
%     % back if necessary
%     [x, nshifts] = shiftdim(x);
%     
%     % check if x is a scalar, and return a row vector if so
%     if isscalar(x)
%       nshifts = 1; 
%     end 
%     
% elseif nargin==3,
%     
%     % if the dimension to take the fft over is given, rearrange the matrix
%     % so this dimension is the first dimension
%     perm = [dim:max(length(size(x)),dim) 1:dim-1];
%     x = permute(x,perm);
% end

% rearrange the matrix so the fft is taken over the first dimension
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = permute(x,perm);

siz = size(x);
[m, n] = size(x);
if ~isscalar(ny) 
  error('MATLAB:interpft:NonScalarN', 'N must be a scalar.'); 
end

% if necessary, increase ny by an integer multiple to make ny > m
% here m is the length of the dimension over which the fft is being taken,
% and ny is the number of points to use including zero padding
if ny > m
   incr = 1;
else
    
    % if ny is zero, return nothing
    if ny == 0
       y = []; 
       return
    end
    
    % otherwise make ny bigger to allow downsampling from integer points
    incr = floor(m/ny) + 1;
    ny = incr*ny;
    
end

% take the fft of the input function
a = fft(x,[],1);

% calculate the nyquist frequency
nyqst = ceil((m+1)/2);

% window the Fourier coefficient if win is given
if win 
    
    % create window using getWin
    wind = getWin(siz(1), win);
    
    % shift the window
    wind = circshift(wind, nyqst);
    
    % create repmat variable
    siz_rep = siz;
    siz_rep(1) = 1;
    
    % repeat
    wind = repmat(wind, siz_rep);

    % apply the window
    a = a.*wind;  
    
end
   
% zero pad with the zeros in the middle
b = [a(1:nyqst,:) ; zeros(ny-m, prod(siz(2:end))) ; a(nyqst+1:m,:)];

if rem(m, 2) == 0 
    
    % if the sequence has an even number of points, make the sequence
    % symmetric, e.g., turn
    % EP P P N   0 0 0 0 0 0   P P
    % into
    % EP P P N/2 0 0 0 0 0 N/2 P P 
    % where EP = end point P = points, N = nyquist point
    b(nyqst,:) = b(nyqst,:)/2;
    b(nyqst+ny-m,:) = b(nyqst,:);
    
end

% take the inverse FFT
y = ifft(b,[],1);

% if the input was real, throw away any residual complex bits
if isreal(x)
    y = real(y); 
end

% make sure the amplitudes are correct by scaling by the new length
y = y * ny / m;

% this gets a downsampled version from an upsampled version
y = y(1:incr:ny,:);  % Skip over extra points when oldny <= m.

% reshape
[y_length, num_signals] = size(y);
y = reshape(y, [y_length, siz(2:end)]);

% sort out the resizing
y = ipermute(y,perm);

% if nargin==2,
%   y = reshape(y,[ones(1,nshifts) size(y,1) siz(2:end)]);
% elseif nargin==3,
%   y = ipermute(y,perm);
% end
