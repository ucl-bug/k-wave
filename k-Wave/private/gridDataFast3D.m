function [tri, bc] = gridDataFast3D(x, y, z, xi, yi, zi)
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

% enforce x, y, z, xi, yi, zi to be column vectors
x = x(:);
y = y(:);
z = z(:);
xi = xi(:);
yi = yi(:);
zi = zi(:);

% check if the MATLAB version is new enough to have delaunayTriangulation
% or DelaunayTri 
if exist('delaunayTriangulation') %#ok<EXIST>
    
    % triangulize the data
    tri = delaunayTriangulation(x, y, z);  
    
elseif exist('DelaunayTri') %#ok<EXIST>
    
    % triangulize the data
    tri = DelaunayTri(x, y, z); %#ok<DDELTRI>
    
else
    
    % give an error
    error('''CartInterp'' set to ''linear'' for 3D simulations is not supported in your version of MATLAB. Try running the simulation with the optional input ''CartInterp'', ''nearest'' or using a different sensor type.');

end

% catch trinagulation error
if isempty(tri)
    error('Data cannot be triangulated.');
end

% find the nearest triangle and the corresponding Barycentric coordinates
[t, bc] = pointLocation(tri, [xi, yi, zi]);

% check points are valid
if any(isnan(t))
    error('Cartesian points must lie within the k-space grid defined by kgrid.');
end

tri = tri(t,:);