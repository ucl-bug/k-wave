% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to extract
%     and display GPU and CPU memory usage. 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 16th July 2013
%     last update - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2017 Bradley Treeby

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

% display current matlab memory usage
if num_outputs == 2 && strncmp(computer, 'PCWIN', 5)
    [mem_usage.user, mem_usage.sys] = memory;
    disp(['  memory used: ' num2str(mem_usage.user.MemUsedMATLAB ./ 1024^3) ' GB (of ' num2str(mem_usage.sys.PhysicalMemory.Total ./ 1024^3) ' GB)']); 
end        

% gpu memory counter for Parallel Computing toolbox
if strcmp(data_cast, 'gpuArray')
    gpu_info = gpuDevice;
    disp(['  GPU memory used: ' num2str((gpu_info.TotalMemory - gpu_info.FreeMemory) ./ 1024^3) ' GB (of ' num2str(gpu_info.TotalMemory ./ 1024^3) ' GB)']);
    mem_usage.gpu.total = gpu_info.TotalMemory;
    mem_usage.gpu.used = gpu_info.TotalMemory - gpu_info.FreeMemory;            
end