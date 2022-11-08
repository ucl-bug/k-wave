function [tri, bc] = gridDataFast2D(x, y, xi, yi)
% Copyright (c) 2012, Chao Huang
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

% Updated 22-02-2017 to use DelaunayTriangulation (Bradley Treeby)

% enforce x, y, xi, yi to be column vectors
x = x(:);
y = y(:);
xi = xi(:);
yi = yi(:);

% check if the MATLAB version is new enough to have delaunayTriangulation
% or DelaunayTri 
if exist('delaunayTriangulation') %#ok<EXIST>
    
    % triangulize the data
    tri = delaunayTriangulation(x, y);  
    
elseif exist('DelaunayTri') %#ok<EXIST>
    
    % triangulize the data
    tri = DelaunayTri(x, y); %#ok<DDELTRI>
    
else
    
    % call the old version of gridDataFast using tsearch
    [tri, bc] = gridDataFast2D_tsearch(x, y, xi, yi);
    return
    
end

% catch trinagulation error
if isempty(tri)
    error('Data cannot be triangulated.');
end

% find the nearest triangle and the corresponding Barycentric coordinates
[t, bc] = pointLocation(tri,[xi yi]);

% check points are valid
if any(isnan(t))
    error('Cartesian points must lie within the k-space grid defined by kgrid.');
end

tri = tri(t,:);