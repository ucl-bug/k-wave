function test_pass = kspaceFirstOrder_compare_1D_2D_3D_plane_waves(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare a plane wave in 1, 2, AS, and 3D to catch any
%     coding bugs between the difference dimensions. 
%
%     84 tests are performed:
%
%         1.  linear + lossless + source.p0 + homogeneous
%         2.  linear + lossless + source.p0 + heterogeneous
%         3.  linear + lossless + source.p (additive) + homogeneous
%         4.  linear + lossless + source.p (additive) + heterogeneous
%         5.  linear + lossless + source.p (additive-no-correction) + homogeneous
%         6.  linear + lossless + source.p (additive-no-correction) + heterogeneous
%         7.  linear + lossless + source.p (dirichlet) + homogeneous
%         8.  linear + lossless + source.p (dirichlet) + heterogeneous
%         9.  linear + lossless + source.u (additive) + homogeneous
%         10. linear + lossless + source.u (additive) + heterogeneous
%         11. linear + lossless + source.u (additive-no-correction) + homogeneous
%         12. linear + lossless + source.u (additive-no-correction) + heterogeneous
%         13. linear + lossless + source.u (dirichlet) + homogeneous
%         14. linear + lossless + source.u (dirichlet) + heterogeneous
%
%         15. linear + lossy + source.p0 + homogeneous
%         16. linear + lossy + source.p0 + heterogeneous
%         17. linear + lossy + source.p (additive) + homogeneous
%         18. linear + lossy + source.p (additive) + heterogeneous
%         19. linear + lossy + source.p (additive-no-correction) + homogeneous
%         20. linear + lossy + source.p (additive-no-correction) + heterogeneous
%         21. linear + lossy + source.p (dirichlet) + homogeneous
%         22. linear + lossy + source.p (dirichlet) + heterogeneous
%         23. linear + lossy + source.u (additive) + homogeneous
%         24. linear + lossy + source.u (additive) + heterogeneous
%         25. linear + lossy + source.u (additive-no-correction) + homogeneous
%         26. linear + lossy + source.u (additive-no-correction) + heterogeneous
%         27. linear + lossy + source.u (dirichlet) + homogeneous
%         28. linear + lossy + source.u (dirichlet) + heterogeneous
%
%         29. linear + stokes + source.p0 + homogeneous
%         30. linear + stokes + source.p0 + heterogeneous
%         31. linear + stokes + source.p (additive) + homogeneous
%         32. linear + stokes + source.p (additive) + heterogeneous
%         33. linear + stokes + source.p (additive-no-correction) + homogeneous
%         34. linear + stokes + source.p (additive-no-correction) + heterogeneous
%         35. linear + stokes + source.p (dirichlet) + homogeneous
%         36. linear + stokes + source.p (dirichlet) + heterogeneous
%         37. linear + stokes + source.u (additive) + homogeneous
%         38. linear + stokes + source.u (additive) + heterogeneous
%         39. linear + stokes + source.u (additive-no-correction) + homogeneous
%         40. linear + stokes + source.u (additive-no-correction) + heterogeneous
%         41. linear + stokes + source.u (dirichlet) + homogeneous
%         42. linear + stokes + source.u (dirichlet) + heterogeneous
%
%         43. nonlinear + lossless + source.p0 + homogeneous
%         44. nonlinear + lossless + source.p0 + heterogeneous
%         45. nonlinear + lossless + source.p (additive) + homogeneous
%         46. nonlinear + lossless + source.p (additive) + heterogeneous
%         47. nonlinear + lossless + source.p (additive-no-correction) + homogeneous
%         48. nonlinear + lossless + source.p (additive-no-correction) + heterogeneous
%         49. nonlinear + lossless + source.p (dirichlet) + homogeneous
%         50. nonlinear + lossless + source.p (dirichlet) + heterogeneous
%         51. nonlinear + lossless + source.u (additive) + homogeneous
%         52. nonlinear + lossless + source.u (additive) + heterogeneous
%         53. nonlinear + lossless + source.u (additive-no-correction) + homogeneous
%         54. nonlinear + lossless + source.u (additive-no-correction) + heterogeneous
%         55. nonlinear + lossless + source.u (dirichlet) + homogeneous
%         56. nonlinear + lossless + source.u (dirichlet) + heterogeneous
%
%         57. nonlinear + lossy + source.p0 + homogeneous
%         58. nonlinear + lossy + source.p0 + heterogeneous
%         59. nonlinear + lossy + source.p (additive) + homogeneous
%         60. nonlinear + lossy + source.p (additive) + heterogeneous
%         61. nonlinear + lossy + source.p (additive-no-correction) + homogeneous
%         62. nonlinear + lossy + source.p (additive-no-correction) + heterogeneous
%         63. nonlinear + lossy + source.p (dirichlet) + homogeneous
%         64. nonlinear + lossy + source.p (dirichlet) + heterogeneous
%         65. nonlinear + lossy + source.u (additive) + homogeneous
%         66. nonlinear + lossy + source.u (additive) + heterogeneous
%         67. nonlinear + lossy + source.u (additive-no-correction) + homogeneous
%         68. nonlinear + lossy + source.u (additive-no-correction) + heterogeneous
%         69. nonlinear + lossy + source.u (dirichlet) + homogeneous
%         70. nonlinear + lossy + source.u (dirichlet) + heterogeneous
%
%         71. nonlinear + stokes + source.p0 + homogeneous
%         72. nonlinear + stokes + source.p0 + heterogeneous
%         73. nonlinear + stokes + source.p (additive) + homogeneous
%         74. nonlinear + stokes + source.p (additive) + heterogeneous
%         75. nonlinear + stokes + source.p (additive-no-correction) + homogeneous
%         76. nonlinear + stokes + source.p (additive-no-correction) + heterogeneous
%         77. nonlinear + stokes + source.p (dirichlet) + homogeneous
%         78. nonlinear + stokes + source.p (dirichlet) + heterogeneous
%         79. nonlinear + stokes + source.u (additive) + homogeneous
%         80. nonlinear + stokes + source.u (additive) + heterogeneous
%         81. nonlinear + stokes + source.u (additive-no-correction) + homogeneous
%         82. nonlinear + stokes + source.u (additive-no-correction) + heterogeneous
%         83. nonlinear + stokes + source.u (dirichlet) + homogeneous
%         84. nonlinear + stokes + source.u (dirichlet) + heterogeneous
%
%     For each test, the plane wave in 2D and 3D is propagated in all
%     possible directions.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 27th July 2011
%     last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011- Bradley Treeby

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

%#ok<*NOPRT>
%#ok<*UNRCH>

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set additional literals to give further permutations of the test
STAGGERED_GRID      = true;
USE_KSPACE          = true;
COMPARISON_THRESH   = 1e-14;
USE_PML             = false;
SMOOTH_P0_SOURCE    = false;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 64;
Nj = 32;
dx = 1;

% define PML properties
PML_size    = 10;
if USE_PML
    PML_alpha_def   = 2;
else
    PML_alpha_def   = 0;
end
    
% define medium properties
c0      = 1500;
rho0    = 1000;
BonA0   = 10;
alpha0  = 5;
y       = 1.5;

% define medium properties
c1      = 2000;
rho1    = 1200;
BonA1   = 5;
alpha1  = 2;

% position of the heterogeneous interface
interface_position = Nx / 2;

% define time array
cfl = 0.1;
Nt  = 600;
dt  = cfl * dx / c0;
t_array = 0:dt:(Nt - 1) * dt;

% define optional inputs
input_args = {'PMLSize', PML_size, 'Smooth', false, 'PlotScale', 'auto', ...
    'UsekSpace', USE_KSPACE, 'UseSG', STAGGERED_GRID, ...
    'PlotSim', plot_simulations};

% define source properties
source_strength = 5e6;
source_position = PML_size - 6;
source_freq     = 200;
source_signal   = source_strength * sin(2 * pi * source_freq * t_array);

% define sensor properties
sensor_position = Nx - PML_size; 

% set pass variable
test_pass = true;

% test names
test_names = {...
    'linear + lossless + source.p0 + homogeneous', ...
    'linear + lossless + source.p0 + heterogeneous', ...
    'linear + lossless + source.p (additive) + homogeneous', ...
    'linear + lossless + source.p (additive) + heterogeneous', ...
    'linear + lossless + source.p (additive-no-correction) + homogeneous', ...
    'linear + lossless + source.p (additive-no-correction) + heterogeneous', ...
    'linear + lossless + source.p (dirichlet) + homogeneous', ...
    'linear + lossless + source.p (dirichlet) + heterogeneous', ...    
    'linear + lossless + source.u (additive) + homogeneous', ...
    'linear + lossless + source.u (additive) + heterogeneous', ...   
    'linear + lossless + source.u (additive-no-correction) + homogeneous', ...
    'linear + lossless + source.u (additive-no-correction) + heterogeneous', ...       
    'linear + lossless + source.u (dirichlet) + homogeneous', ...
    'linear + lossless + source.u (dirichlet) + heterogeneous', ...        
    'linear + lossy + source.p0 + homogeneous', ...
    'linear + lossy + source.p0 + heterogeneous', ...
    'linear + lossy + source.p (additive) + homogeneous', ...
    'linear + lossy + source.p (additive) + heterogeneous', ...
    'linear + lossy + source.p (additive-no-correction) + homogeneous', ...
    'linear + lossy + source.p (additive-no-correction) + heterogeneous', ...    
    'linear + lossy + source.p (dirichlet) + homogeneous', ...
    'linear + lossy + source.p (dirichlet) + heterogeneous', ...    
    'linear + lossy + source.u (additive) + homogeneous', ...
    'linear + lossy + source.u (additive) + heterogeneous', ...   
    'linear + lossy + source.u (additive-no-correction) + homogeneous', ...
    'linear + lossy + source.u (additive-no-correction) + heterogeneous', ...       
    'linear + lossy + source.u (dirichlet) + homogeneous', ...
    'linear + lossy + source.u (dirichlet) + heterogeneous', ...     
    'linear + stokes + source.p0 + homogeneous', ...
    'linear + stokes + source.p0 + heterogeneous', ...
    'linear + stokes + source.p (additive) + homogeneous', ...
    'linear + stokes + source.p (additive) + heterogeneous', ...
    'linear + stokes + source.p (additive-no-correction) + homogeneous', ...
    'linear + stokes + source.p (additive-no-correction) + heterogeneous', ...    
    'linear + stokes + source.p (dirichlet) + homogeneous', ...
    'linear + stokes + source.p (dirichlet) + heterogeneous', ...    
    'linear + stokes + source.u (additive) + homogeneous', ...
    'linear + stokes + source.u (additive) + heterogeneous', ...    
    'linear + stokes + source.u (additive-no-correction) + homogeneous', ...
    'linear + stokes + source.u (additive-no-correction) + heterogeneous', ...        
    'linear + stokes + source.u (dirichlet) + homogeneous', ...
    'linear + stokes + source.u (dirichlet) + heterogeneous', ...      
    'nonlinear + lossless + source.p0 + homogeneous', ...
    'nonlinear + lossless + source.p0 + heterogeneous', ...
    'nonlinear + lossless + source.p (additive) + homogeneous', ...
    'nonlinear + lossless + source.p (additive) + heterogeneous', ...
    'nonlinear + lossless + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + lossless + source.p (additive-no-correction) + heterogeneous', ...    
    'nonlinear + lossless + source.p (dirichlet) + homogeneous', ...
    'nonlinear + lossless + source.p (dirichlet) + heterogeneous', ...    
    'nonlinear + lossless + source.u (additive) + homogeneous', ...
    'nonlinear + lossless + source.u (additive) + heterogeneous', ...   
    'nonlinear + lossless + source.u (additive-no-correction) + homogeneous', ...
    'nonlinear + lossless + source.u (additive-no-correction) + heterogeneous', ...      
    'nonlinear + lossless + source.u (dirichlet) + homogeneous', ...
    'nonlinear + lossless + source.u (dirichlet) + heterogeneous', ...       
    'nonlinear + lossy + source.p0 + homogeneous', ...
    'nonlinear + lossy + source.p0 + heterogeneous', ...
    'nonlinear + lossy + source.p (additive) + homogeneous', ...
    'nonlinear + lossy + source.p (additive) + heterogeneous', ...
    'nonlinear + lossy + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + lossy + source.p (additive-no-correction) + heterogeneous', ...    
    'nonlinear + lossy + source.p (dirichlet) + homogeneous', ...
    'nonlinear + lossy + source.p (dirichlet) + heterogeneous', ...    
    'nonlinear + lossy + source.u (additive) + homogeneous', ...
    'nonlinear + lossy + source.u (additive) + heterogeneous', ...    
    'nonlinear + lossy + source.u (additive-no-correction) + homogeneous', ...
    'nonlinear + lossy + source.u (additive-no-correction) + heterogeneous', ...      
    'nonlinear + lossy + source.u (dirichlet) + homogeneous', ...
    'nonlinear + lossy + source.u (dirichlet) + heterogeneous', ... 
    'nonlinear + stokes + source.p0 + homogeneous', ...
    'nonlinear + stokes + source.p0 + heterogeneous', ...
    'nonlinear + stokes + source.p (additive) + homogeneous', ...
    'nonlinear + stokes + source.p (additive) + heterogeneous', ...
    'nonlinear + stokes + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + stokes + source.p (additive-no-correction) + heterogeneous', ...    
    'nonlinear + stokes + source.p (dirichlet) + homogeneous', ...
    'nonlinear + stokes + source.p (dirichlet) + heterogeneous', ...    
    'nonlinear + stokes + source.u (additive) + homogeneous', ...
    'nonlinear + stokes + source.u (additive) + heterogeneous', ...    
    'nonlinear + stokes + source.u (additive-no-correction) + homogeneous', ...
    'nonlinear + stokes + source.u (additive-no-correction) + heterogeneous', ...     
    'nonlinear + stokes + source.u (dirichlet) + homogeneous', ...
    'nonlinear + stokes + source.u (dirichlet) + heterogeneous', ...     
    };

% lists used to set properties
p0_tests = [1, 2, 15, 16, 29, 30, 43, 44, 57, 58, 71, 72];
p_tests  = [3:8, 17:22, 31:36, 45:50, 59:64, 73:78];
u_tests  = [9:14, 23:28, 37:42, 51:56, 65:70, 79:84];
additive_tests = [3,  4,  9,  10, 17, 18, 23, 24, 31, 32, 37, 38, ...
                  45, 46, 51, 52, 59, 60, 65, 66, 73, 74, 79, 80];
additive_no_correction_tests = additive_tests + 2;
dirichlet_tests = additive_tests + 4;

% =========================================================================
% SIMULATIONS
% =========================================================================

% loop through tests
for test_num = 1:84

    % clear structures
    clear source medium
    
    % update command line
    disp(['Running Test: ' test_names{test_num}]);
       
    % assign medium properties 
    medium.sound_speed  = c0;
    medium.density      = rho0;
    PML_alpha = PML_alpha_def;
    switch test_num
        case {1,2,3,4,5,6,7,8,9,10,11,12,13,14}
            
            % linear + lossless
            
        case {15,16,17,18,19,20,21,22,23,24,25,26,27,28}
            
            % linear + lossy
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = y;    
            
        case {29,30,31,32,33,34,35,36,37,38,39,40,41,42}
            
            % linear + stokes
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = 2;    
            medium.alpha_mode   = 'stokes';
            
        case {43,44,45,46,47,48,49,50,51,52,53,54,55,56}
            
            % nonlinear + lossless
            medium.BonA         = BonA0;
            
        case {57,58,59,60,61,62,63,64,65,66,67,68,69,70}
            
            % nonlinear + lossy
            medium.BonA         = BonA0;
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = y;
            
        case {71,72,73,74,75,76,77,78,79,80,81,82,83,84}
            
            % nonlinear + stokes
            medium.BonA         = BonA0;
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = 2;
            medium.alpha_mode   = 'stokes';
            
        otherwise
            
            error('Unknown test number.');
            
    end
    
    % ----------------
    % 1D SIMULATION: X
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nx, dx);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, 1, 1, 1)
    end
    
    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, 1);
        source.p0(source_position) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true); 
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nx, 1);
        source.p_mask(source_position) = 1;
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end     
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, 1);
        source.u_mask(source_position) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end     
    else
        error('Unknown source condition.');
    end
    
    % sensor
    sensor.mask = zeros(Nx, 1);
    sensor.mask(sensor_position) = 1;

    % run simulation
    sensor_data_1D = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', PML_alpha);
    
    % ----------------
    % 2D SIMULATION: X
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Nj, dx);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Nj, 1, 1)
    end    

    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Nj);
        source.p0(source_position, :) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nx, Nj);
        source.p_mask(source_position, :) = 1; 
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end     
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Nj);
        source.u_mask(source_position, :) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end      
    else
        error('Unknown source condition.');
    end
    
    % sensor
    sensor.mask = zeros(Nx, Nj);
    sensor.mask(sensor_position, Nj/2) = 1;

    % run simulation
    sensor_data_2D_x = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', [PML_alpha, 0]);
    
    % ----------------
    % 2D SIMULATION: Y
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nj, dx, Nx, dx);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nj, Nx, 1, 2)
    end      
    
    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nj, Nx);
        source.p0(:, source_position) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nj, Nx);
        source.p_mask(:, source_position) = 1; 
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end    
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nj, Nx);
        source.u_mask(:, source_position) = 1;
        source.uy = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end     
    else
        error('Unknown source condition.');
    end
    
    % sensor
    sensor.mask = zeros(Nj, Nx);
    sensor.mask(Nj/2, sensor_position) = 1;

    % run simulation
    sensor_data_2D_y = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', [0, PML_alpha]);    

    % -------------------------------------
    % 2D AXISYMMETRIC SIMULATION: X
    % -------------------------------------

    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Nj, dx);
    kgrid.t_array = t_array;
    
    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Nj, 1, 1)
    end    

    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Nj);
        source.p0(source_position, :) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nx, Nj);
        source.p_mask(source_position, :) = 1; 
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end      
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Nj);
        source.u_mask(source_position, :) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end   
    else
        error('Unknown source condition.');
    end
    
    % sensor
    sensor.mask = zeros(Nx, Nj);
    sensor.mask(sensor_position, Nj/2) = 1;

    % run simulation
    sensor_data_AS_wsws_x = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:}, ...
        'RadialSymmetry', 'WSWS-FFT', 'PMLAlpha', [PML_alpha, 0]);   
    
    % ----------------
    % 3D SIMULATION: X
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Nj, dx, Nj, dx);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nx, Nj, Nj, 1)
    end      
    
    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Nj, Nj);
        source.p0(source_position, :, :) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nx, Nj, Nj);
        source.p_mask(source_position, :, :) = 1;
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end    
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Nj, Nj);
        source.u_mask(source_position, :, :) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end         
    else
        error('Unknown source condition.');
    end

    % sensor
    sensor.mask = zeros(Nx, Nj, Nj);
    sensor.mask(sensor_position, Nj/2, Nj/2) = 1;

    % run simulation
    sensor_data_3D_x = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', [PML_alpha, 0, 0]);
    
    % ----------------
    % 3D SIMULATION: Y
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nj, dx, Nx, dx, Nj, dx);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nj, Nx, Nj, 2)
    end          
    
    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nj, Nx, Nj);
        source.p0(:, source_position, :) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nj, Nx, Nj);
        source.p_mask(:, source_position, :) = 1; 
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end        
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nj, Nx, Nj);
        source.u_mask(:, source_position, :) = 1;
        source.uy = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end      
    else
        error('Unknown source condition.');        
    end

    % sensor
    sensor.mask = zeros(Nj, Nx, Nj);
    sensor.mask(Nj/2, sensor_position, Nj/2) = 1;

    % run simulation
    sensor_data_3D_y = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', [0, PML_alpha, 0]);
    
    % ----------------
    % 3D SIMULATION: Z
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nj, dx, Nj, dx, Nx, dx);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties(Nj, Nj, Nx, 3)
    end          
    
    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nj, Nj, Nx);
        source.p0(:, :, source_position) = source_strength;
        if SMOOTH_P0_SOURCE
            source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nj, Nj, Nx);
        source.p_mask(:, :, source_position) = 1;
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nj, Nj, Nx);
        source.u_mask(:, :, source_position) = 1;
        source.uz = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end            
    else
        error('Unknown source condition.');        
    end

    % sensor
    sensor.mask = zeros(Nj, Nj, Nx);
    sensor.mask(Nj/2, Nj/2, sensor_position) = 1;

    % run simulation
    sensor_data_3D_z = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', [0, 0, PML_alpha]);
    
    % -------------
    % COMPARISON
    % -------------

    % compute errors
    ref_max = max(abs(sensor_data_1D(:)));
    diff_1D_2D_x        = max(abs(sensor_data_1D(:) - sensor_data_2D_x(:))) / ref_max
    diff_1D_2D_y        = max(abs(sensor_data_1D(:) - sensor_data_2D_y(:))) / ref_max
    diff_1D_AS_wsws_x   = max(abs(sensor_data_1D(:) - sensor_data_AS_wsws_x(:))) / ref_max
    diff_1D_3D_x        = max(abs(sensor_data_1D(:) - sensor_data_3D_x(:))) / ref_max
    diff_1D_3D_y        = max(abs(sensor_data_1D(:) - sensor_data_3D_y(:))) / ref_max
    diff_1D_3D_z        = max(abs(sensor_data_1D(:) - sensor_data_3D_z(:))) / ref_max
    
    % ignore error from axisymemtric code for lossy case
    if isfield(medium, 'alpha_power') && (medium.alpha_power == y)
        diff_1D_AS_wsws_x = 0;
    end
       
    % check for test pass
    if (diff_1D_2D_x > COMPARISON_THRESH) || ...
       (diff_1D_2D_y > COMPARISON_THRESH) || ...
       (diff_1D_AS_wsws_x > COMPARISON_THRESH) || ...
       (diff_1D_3D_x > COMPARISON_THRESH) || ...
       (diff_1D_3D_y > COMPARISON_THRESH) || ...
       (diff_1D_3D_z > COMPARISON_THRESH)
        test_pass = false;
    end
    
    % -------------
    % PLOTTING
    % -------------    
    
    if plot_comparisons
        figure;
        
        subplot(4, 2, 1), plot(sensor_data_1D);
        title({test_names{test_num}; '1D'});
        
        subplot(4, 2, 2), plot(sensor_data_2D_x);
        title(['2D X, L_{inf} = ' num2str(diff_1D_2D_x)]);
        
        subplot(4, 2, 3), plot(sensor_data_2D_y);
        title(['2D Y, L_{inf} = ' num2str(diff_1D_2D_y)]);
        
        subplot(4, 2, 4), plot(sensor_data_AS_wsws_x);
        title(['AS WSWS, L_{inf} = ' num2str(diff_1D_AS_wsws_x)]);  
              
        subplot(4, 2, 5), plot(sensor_data_3D_x);
        title(['3D X, L_{inf} = ' num2str(diff_1D_3D_x)]);
        
        subplot(4, 2, 6), plot(sensor_data_3D_y);
        title(['3D Y, L_{inf} = ' num2str(diff_1D_3D_y)]);
        
        subplot(4, 2, 7), plot(sensor_data_3D_z);
        title(['3D Z, L_{inf} = ' num2str(diff_1D_3D_z)]);
        
        scaleFig(1, 1.5);
        drawnow;
    end
    
end

% =========================================================================
% SUB FUNCTIONS
% =========================================================================

function setMaterialProperties(N1, N2, N3, direction)

    % sound speed and density
    medium.sound_speed = c0 * ones(N1, N2, N3);
    medium.density = rho0 * ones(N1, N2, N3);
    switch direction
        case 1
            medium.sound_speed(interface_position:end, :, :) = c1;
            medium.density(interface_position:end, :, :) = rho1;
        case 2
            medium.sound_speed(:, interface_position:end, :) = c1;
            medium.density(:, interface_position:end, :) = rho1;
        case 3
            medium.sound_speed(:, :, interface_position:end) = c1;
            medium.density(:, :, interface_position:end) = rho1;
    end

    % absorption
    if isfield(medium, 'alpha_coeff')
        medium.alpha_coeff = alpha0 * ones(N1, N2, N3);
        switch direction
            case 1
                medium.alpha_coeff(interface_position:end, :, :) = alpha1;
            case 2
                medium.alpha_coeff(:, interface_position:end, :) = alpha1;
            case 3
                medium.alpha_coeff(:, :, interface_position:end) = alpha1;
        end

    end

    % nonlinearity
    if isfield(medium, 'BonA')
        medium.BonA = BonA0 * ones(N1, N2, N3);
        switch direction
            case 1
                medium.BonA(interface_position:end, :, :) = BonA1;
            case 2
                medium.BonA(:, interface_position:end, :) = BonA1;
            case 3
                medium.BonA(:, :, interface_position:end) = BonA1;
        end

    end        
    
end

end