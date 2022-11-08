function line = makeLine(Nx, Ny, startpoint, varargin)
%MAKELINE Create a binary map of a straight line within a 2D grid.
%
% DESCRIPTION:
%     makeLine creates a binary map of a line within a two-dimensional grid
%     (the line is denoted by 1's in the matrix with 0's elsewhere). The
%     line can be defined either by its endpoints or by a starting point,
%     an angle (direction), and a length.
%
% USAGE:
%     line = makeLine(Nx, Ny, startpoint, endpoint)
%     line = makeLine(Nx, Ny, startpoint, angle, length)
%
% INPUTS:
%     Nx, Ny          - size of the 2D grid [grid points]
%     startpoint      - two-element vector identifying the start point
%                       for the line [grid points]  
%     endpoint        - two-element vector identifying the end
%                       point for the line [grid points] 
%     angle           - the direction of the line [radians]
%     length          - length of the line [grid points]
%
% OUTPUTS:
%     line            - 2D binary map of a straight line
%
% ABOUT:
%     author          - Ben Cox
%     date            - 4th September 2012
%     last update     - 13th January 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2019 Ben Cox
%
% See also makeArc, makeCircle, makeDisc

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

% startpoint must be a two-element row vector identifying a point within
% the grid

% =========================================================================
% INPUT CHECKING
% =========================================================================

if length(startpoint) ~= 2
    error('startpoint should be a two-element vector.')
end
if size(startpoint, 2) ~= 2
    startpoint = startpoint.';
end
if (sum(startpoint < 1) > 0) || (startpoint(1) > Nx) || (startpoint(2) > Ny)
   error('The starting point must lie within the grid, between [1 1] and [Nx Ny].') 
end

% =========================================================================
% LINE BETWEEN TWO POINTS OR ANGLED LINE?
% =========================================================================

switch size(varargin, 2)
    case 1 
        linetype = 'AtoB';
        a = startpoint;
        b = varargin{1};
    case 2
        linetype = 'angled';
        angle = varargin{1};
        linelength = varargin{2};
    otherwise
        error('There must be four or five inputs.')
end

% =========================================================================
% MORE INPUT CHECKING
% =========================================================================

if strcmp(linetype, 'AtoB')
    
    % a and b must be different points
    if (a(1) == b(1)) && (a(2) == b(2))
        error('The first and last points cannot be the same.')
    end

    % end point must be a two-element row vector 
    if length(b) ~= 2
        error('endpoint should be a two-element vector.')
    end
    if size(b, 2) ~= 2
        b = b.';
    end
    
    % a and b must be within the grid
    xx = [a(1), b(1)];
    yy = [a(2), b(2)];
    if (sum(a < 1) > 0) || (sum(b < 1) > 0) || (sum(xx > Nx) > 0) || (sum(yy > Ny) > 0)
       error('Both the start and end points must lie within the grid.') 
    end

    % a and b must be integers
    if sum(rem([a, b], 1)) > 0
        error('Both the start and end points must be given as pairs of integers.')
    end
    
end
    
if strcmp(linetype, 'angled')
    
    % angle must lie between -pi and pi 
    angle = rem(angle, 2*pi);
    if angle > pi
        angle = angle - 2*pi;
    elseif angle < -pi
        angle = angle + 2*pi;
    end
    
end


% =========================================================================
% CALCULATE A LINE FROM A TO B
% =========================================================================

if strcmp(linetype, 'AtoB')
    
    % define an empty grid to hold the line
    line = zeros(Nx, Ny);

    % find the equation of the line
    m = (b(2) - a(2)) / (b(1) - a(1));    % gradient of the line
    c = a(2) - m * a(1);                  % where the line crosses the y axis

    if abs(m) < 1

        % start at the end with the smallest value of x
        if a(1) < b(1)
            x = a(1);
            y = a(2);
            x_end = b(1);
        else
            x = b(1);
            y = b(2);
            x_end = a(1);
        end

        % fill in the first point
        line(x, y) = 1;

        while x<x_end

            % next points to try are
            poss_x = [x   x   x+1   x+1 x+1];
            poss_y = [y-1 y+1 y-1   y   y+1];

            % find the point closest to the line
            true_y = m * poss_x + c;
            diff = (poss_y - true_y).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));
            
            % add the point to the line
            line(x, y) = 1;

        end

    elseif ~isinf(abs(m))

        % start at the end with the smallest value of y
        if a(2) < b(2)
            x = a(1);
            y = a(2);
            y_end = b(2);
        else
            x = b(1);
            y = b(2);
            y_end = a(2);
        end

        % fill in the first point
        line(x, y) = 1;

        while y < y_end

            % next points to try are
            poss_y = [y   y   y+1   y+1 y+1];
            poss_x = [x-1 x+1 x-1   x   x+1];

            % find the point closest to the line
            true_x = (poss_y - c) / m;
            diff = (poss_x - true_x).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));

            % add the point to the line
            line(x, y) = 1;

        end

    else % m = +-Inf

        % start at the end with the smallest value of y
        if a(2) < b(2)
            x = a(1);
            y = a(2);
            y_end = b(2);
        else
            x = b(1);
            y = b(2);
            y_end = a(2);
        end

        % fill in the first point
        line(x, y) = 1;

        while y < y_end

            % next point
            y = y + 1;

            % add the point to the line
            line(x, y) = 1;

        end   

    end
        
% =========================================================================
% CALCULATE AN ANGLED LINE
% =========================================================================

elseif strcmp(linetype, 'angled')
    
    % define an empty grid to hold the line
    line = zeros(Nx, Ny);

    % start at the atart
    x = startpoint(1);
    y = startpoint(2);

    % fill in the first point
    line(x, y) = 1;

    % initialise the current length of the line
    line_length = 0;

    if abs(angle) == pi

        while line_length < linelength

            % next point
            y = y + 1;

            % stop the points incrementing at the edges
            if y > Ny
                break;
            end
                
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end   

    elseif (angle < pi) && (angle > pi/2)

        % define the equation of the line
        m = -tan(angle - pi/2);   % gradient of the line
        c = y - m * x;            % where the line crosses the y axis

        while line_length<linelength

            % next points to try are
            poss_x = [x-1 x-1 x  ];
            poss_y = [y   y+1 y+1];

            % find the point closest to the line
            true_y = m * poss_x + c;
            diff = (poss_y - true_y).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));

            % stop the points incrementing at the edges
            if (x < 1) || (y > Ny)
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end

    elseif angle == pi/2

        while line_length < linelength

            % next point
            x = x - 1;

            % stop the points incrementing at the edges
            if x < 1
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end       

    elseif (angle < pi/2) && (angle > 0)

        % define the equation of the line
        m = tan(pi/2 - angle);    % gradient of the line
        c = y - m * x;            % where the line crosses the y axis

        while line_length < linelength

            % next points to try are
            poss_x = [x-1 x-1 x ];
            poss_y = [y   y-1 y-1];

            % find the point closest to the line
            true_y = m * poss_x + c;
            diff = (poss_y - true_y).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));

            % stop the points incrementing at the edges
            if (x < 1) || (y < 1)
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end

    elseif angle==0

        while line_length < linelength

            % next point
            y = y - 1;

            % stop the points incrementing at the edges
            if y < 1
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end       

    elseif (angle < 0) && (angle > -pi/2)

        % define the equation of the line
        m = -tan(pi/2 + angle);   % gradient of the line
        c = y - m * x;            % where the line crosses the y axis

        while line_length < linelength

            % next points to try are
            poss_x = [x+1 x+1 x ];
            poss_y = [y   y-1 y-1];

            % find the point closest to the line
            true_y = m * poss_x + c;
            diff = (poss_y - true_y).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));

            % stop the points incrementing at the edges
            if (x > Nx) || (y < 1)
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end    

    elseif angle == -pi/2

        while line_length < linelength

            % next point
            x = x + 1;

            % stop the points incrementing at the edges
            if x > Nx
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end       

    elseif (angle < -pi/2) && (angle > -pi)

        % define the equation of the line
        m = tan(-angle - pi/2);     % gradient of the line
        c = y - m * x;              % where the line crosses the y axis

        while line_length < linelength

            % next points to try are
            poss_x = [x+1 x+1  x ];
            poss_y = [y   y+1 y+1];

            % find the point closest to the line
            true_y = m * poss_x + c;
            diff = (poss_y - true_y).^2;
            index = find(diff == min(diff));

            % the next point
            x = poss_x(index(1));
            y = poss_y(index(1));

            % stop the points incrementing at the edges
            if (x > Nx) || (y > Ny)
                break;
            end
            
            % add the point to the line
            line(x, y) = 1;

            % calculate the current length of the line
            line_length = sqrt((x - startpoint(1))^2 + (y - startpoint(2))^2);

        end    
    end
end