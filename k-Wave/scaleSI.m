function [x_sc, scale, prefix, prefix_fullname] = scaleSI(x)
%SCALESI Scale a number to nearest SI unit prefix.
%
% DESCRIPTION:
%     scaleSI scales the input x to use the nearest SI unit prefix while
%     keeping 1000 > x > 1. 
%
%     For example, scaleSI(0.00001) would give x_sc = '10u', scale = 1e6,
%     prefix = 'u', and prefix_fullname = 'micro'
%
% USAGE:
%     x_sc = scaleSI(x)
%     [x_sc, scale] = scaleSI(x)
%     [x_sc, scale, prefix] = scaleSI(x)
%     [x_sc, scale, prefix, prefix_fullname] = scaleSI(x)
%
% INPUTS:
%     x               - number to scale
%
% OUTPUTS:
%     x_sc            - string of scaled input and prefix
%     scale           - numeric scale factor
%     prefix          - single character scale prefix
%     prefix_fullname - full SI name for prefix
%       
% ABOUT:
%     author          - Bradley Treeby
%     date            - 15th June 2009
%     last update     - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also scaleTime

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

% force the input to be a scalar
x = max(x(:));

% check for a negative input
if x < 0
    x = -x;
    negative = true;
else
    negative = false;
end

if x == 0
    
    % if x is zero, don't scale
    x_sc = x;
    prefix = '';
    prefix_fullname = '';
    scale = 1;
    
elseif x < 1
    
    % update index and input
    x_sc = x * 1e3;
    sym_index = 1;
       
    % find scaling parameter
    while x_sc < 1 && sym_index < 8
        x_sc = x_sc * 1e3;
        sym_index = sym_index + 1;
    end

    % define SI unit scalings
    switch sym_index
        case 1
            prefix = 'm';
            prefix_fullname = 'milli';
            scale = 1e3;
        case 2
            prefix = 'u';
            prefix_fullname = 'micro';
            scale = 1e6;
        case 3
            prefix = 'n';
            prefix_fullname = 'nano';
            scale = 1e9;
        case 4
            prefix = 'p';
            prefix_fullname = 'pico';
            scale = 1e12;
        case 5
            prefix = 'f';
            prefix_fullname = 'femto';
            scale = 1e15;
        case 6
            prefix = 'a';
            prefix_fullname = 'atto';
            scale = 1e18;
        case 7
            prefix = 'z';
            prefix_fullname = 'zepto';
            scale = 1e21;
        case 8
            prefix = 'y';
            prefix_fullname = 'yocto';
            scale = 1e24;
    end    
    
elseif x >= 1000
    
    % update index and input
    x_sc = x * 1e-3;
    sym_index = 1;
    
    % find scaling parameter
    while x_sc >= 1000 && sym_index < 8
        x_sc = x_sc * 1e-3;
        sym_index = sym_index + 1;
    end
        
    % define SI unit scalings
    switch sym_index
        case 1
            prefix = 'k';
            prefix_fullname = 'kilo';
            scale = 1e-3;
        case 2
            prefix = 'M';
            prefix_fullname = 'mega';
            scale = 1e-6;
        case 3
            prefix = 'G';
            prefix_fullname = 'giga';
            scale = 1e-9;
        case 4
            prefix = 'T';
            prefix_fullname = 'tera';
            scale = 1e-12;
        case 5
            prefix = 'P';
            prefix_fullname = 'peta';
            scale = 1e-15;
        case 6
            prefix = 'E';
            prefix_fullname  = 'exa';
            scale = 1e-18;
        case 7
            prefix = 'Z';
            prefix_fullname = 'zetta';
            scale = 1e-21;
        case 8
            prefix = 'Y';
            prefix_fullname = 'yotta';
            scale = 1e-24;
    end
    
else
    
    % if x is between 1 and 1000, don't scale
    x_sc = x;
    prefix = '';
    prefix_fullname = '';
    scale = 1;
    
end

% form scaling into a string
if negative
    x_sc = ['-', num2str(x_sc), prefix];
else
    x_sc = [num2str(x_sc), prefix];
end