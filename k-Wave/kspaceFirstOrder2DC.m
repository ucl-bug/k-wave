function sensor_data = kspaceFirstOrder2DC(varargin)
%KSPACEFIRSTORDER2DC 2D time-domain simulation of wave propagation using C++ code.
%
% DESCRIPTION:
%     kspaceFirstOrder2DC provides a blind interface to the C++ version of
%     kspaceFirstOrder2D (called kspaceFirstOrder-OMP) in the same way as
%     kspaceFirstOrder3DC. Note, the C++ code does not support all input
%     options, and all display options are ignored (only command line
%     outputs are given). See the k-Wave user manual for more information.
%
%     The function works by appending the optional input 'SaveToDisk' to
%     the user inputs and then calling kspaceFirstOrder2D to save the input
%     files to disk. The contents of sensor.record (if set) are parsed as
%     input flags, and the C++ code is run using the system command. The
%     output files are then automatically loaded from disk and returned in
%     the same fashion as kspaceFirstOrder2D. The input and output files
%     are saved to the temporary directory native to the operating system,
%     and are deleted after the function runs.
%
%     For small simulations, running the simulation on a smaller number of
%     cores can improve performance as the matrices are often small enough
%     to fit within cache. It is recommended to adjust the value of
%     'NumThreads' to optimise performance for a given simulation size and
%     computer hardware. By default, simulations smaller than 128^2 are
%     set to run using a single thread (this behaviour can be over-ridden
%     using the 'NumThreads' option). In some circumstances, for very small
%     simulations, the C++ code can be slower than the MATLAB code.
%
%     This function requires the C++ binary/executable of 
%     kspaceFirstOrder-OMP to be downloaded from
%     http://www.k-wave.org/download.php and placed in the "binaries"
%     directory of the k-Wave toolbox (the same binary is used for
%     simulations in 2D, 3D, and axisymmetric coordinates). Alternatively,
%     the name and location of the binary can be  specified using the
%     optional input parameters 'BinaryName' and 'BinariesPath'.
%
% USAGE:
%     see kspaceFirstOrder2D, kspaceFirstOrder2DG, kspaceFirstOrder3DC
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 28th February 2019
%     last update     - 10th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby
%
% See also kspaceFirstOrder2D, kspaceFirstOrder3DC

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

% This function is essentially a wrapper and directly uses the capabilities
% of kspaceFirstOrder3DC by replacing the binary name with the name of the
% GPU binary. 

% Check for the function name input. If not defined, set the default name
% of the k-Wave MATLAB function to use to generate input file.
if ~any(strcmp('FunctionName', varargin))
    varargin = [varargin {'FunctionName', 'kspaceFirstOrder2D'}];
end

% Check for the number of thread inputs. If not defined and the simulation
% is smaller than 128^2, run the simulation using a single thread. Note,
% this doesn't account for the PML being outside the grid.
if ~any(strcmp('NumThreads', varargin))
    
    % assign grid input
    kgrid = varargin{1};
    
    % check grid size and set threads option
    if (kgrid.Nx * kgrid.Ny) < (128^2)
        varargin = [varargin {'NumThreads', 1}];
    end
    
end

% pass inputs to kspaceFirstOrder3DC
sensor_data = kspaceFirstOrder3DC(varargin{:});