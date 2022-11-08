function [tri, bc] = gridDataFast2D_tsearch(x, y, xi, yi)
% DESCRIPTION:
%     Stripped down version of GRIDDATA that removes inbuilt data checking
%     and allows input and output of the delaunay triangulation for use on
%     subsequent calls with the exact same set of data coordinates.
%
% USAGE:
%     [tri, bc] = gridDataFast(x, y, xi, yi)
%
% INPUT/OUTPUT
%     x, y, xi, yi, are the same as griddata
%
% ABOUT:
%     author      - John Wilkin
%     date        - October 2002
%     modified by - Bradley Treeby
%     last update - 4th October 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2002-2017 John Wilkin and Bradley Treeby

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

% enforce x,y and z to be column vectors
sz = numel(x);
x = reshape(x, sz, 1);
y = reshape(y, sz, 1);
xi = xi(:); 
yi = yi(:);
x = x(:); 
y = y(:);

% triangulize the data
tri = delaunayn([x, y]);

% catch trinagulation error
if isempty(tri)    
    error('Data cannot be triangulated.');
end

% find the nearest triangle (t)
t = tsearch(x, y, tri, xi, yi);

% only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out)
    t(out) = ones(size(out));
end  
tri = tri(t,:);

% compute Barycentric coordinates
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1)))...
    - (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
bc(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi)...
    - (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
bc(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi)...
    - (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
bc(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi)...
    - (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
bc(out,:) = zeros(length(out),3);