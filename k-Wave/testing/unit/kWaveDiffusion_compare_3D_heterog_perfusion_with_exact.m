function test_pass = kWaveDiffusion_compare_3D_heterog_perfusion_with_exact(plot_comparisons, plot_simulations) 
% DESCRIPTION:
%       Unit test to compare the k-Wave diffusion class with an exact
%       solution for Pennes' bioheat equation.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 16th August 2015
%       last update - 7th August 2023
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2023 Bradley Treeby

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

% set pass variable
test_pass = true;

% set comparison threshold
COMPARISON_THRESH = 3e-5;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% run the whole test twice, once with the 'DataCast' option
for test_ind = 1:2

    clear kgrid medium source

    % =====================================================================
    % SIMULATION
    % =====================================================================
    
    % create grid
    Nx = 64;
    Ny = 64;
    Nz = 64;
    dx = 1e-3;
    dy = 1e-3;
    dz = 1e-3;
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
    
    % define medium properties
    density_1              = 1079;  % [kg/m^3]
    thermal_conductivity_1 = 0.52;  % [W/(m.K)]
    specific_heat_1        = 3540;  % [J/(kg.K)]
    
    density_2              = 1200;  % [kg/m^3]
    thermal_conductivity_2 = 0.13;  % [W/(m.K)]
    specific_heat_2        = 4000;  % [J/(kg.K)]
    
    % define blood properties
    blood_density          = 1060;  % [kg/m^3]
    blood_specific_heat    = 3617;  % [J/(kg.K)]
    blood_ambient_temp     = 37;    % [degC]
    blood_perfusion_rate_1 = 0.005;  % [1/s]
    blood_perfusion_rate_2 = 0.02;  % [1/s]
    
    % calculate perfusion
    perfusion_coeff_1 = (blood_density .* blood_perfusion_rate_1 .* blood_specific_heat) ./ (density_1 .* specific_heat_1);
    perfusion_coeff_2 = (blood_density .* blood_perfusion_rate_2 .* blood_specific_heat) ./ (density_2 .* specific_heat_2);
    
    % define heterogeneous medium properties
    medium.density              = density_1 .* ones(Nx, Ny, Nz);
    medium.thermal_conductivity = thermal_conductivity_1 .* ones(Nx, Ny, Nz);
    medium.specific_heat        = specific_heat_1 .* ones(Nx, Ny, Nz);
    medium.perfusion_coeff      = perfusion_coeff_1 .* ones(Nx, Ny, Nz);
    
    medium.density(kgrid.y > 0)              = density_2;
    medium.thermal_conductivity(kgrid.y > 0) = thermal_conductivity_2;
    medium.specific_heat(kgrid.y > 0)        = specific_heat_2;
    medium.perfusion_coeff(kgrid.y > 0)      = perfusion_coeff_2;
    
    % define arterial temperature
    medium.blood_ambient_temperature = blood_ambient_temp;     	% [degC]
    
    % set initial temperature distribution
    s1 = smooth(makeBall(Nx, Ny, Nz, Nx/2, Ny/4, Nz/2, 3));
    s2 = smooth(makeBall(Nx, Ny, Nz, Nx/2, 3*Ny/4, Nz/2, 3));
    
    % define source term
    source.T0  = 37 + s1 + s2;
    
    % set number of time steps
    Nt = 75;
    dt = 0.5;
    
    % set input args
    switch test_ind
        case 1
            input_args = {'PlotSim', plot_simulations, 'PlotScale', [37, 38]};
        case 2
            input_args = {'PlotSim', plot_simulations, 'PlotScale', [37, 38], 'DataCast', 'gpuArray-single'};
    end
    
    % create kWaveDiffusion object and take time steps
    medium.diffusion_coeff_ref = 'max';
    medium.perfusion_coeff_ref = 'min';
    kdiff_max = kWaveDiffusion(kgrid, medium, source, [], input_args{:});
    kdiff_max.takeTimeStep(Nt, dt);
    
    % create kWaveDiffusion object again using opposite reference value and
    % take time steps 
    medium.diffusion_coeff_ref = 'min';
    medium.perfusion_coeff_ref = 'max';
    kdiff_min = kWaveDiffusion(kgrid, medium, source, [], input_args{:});
    kdiff_min.takeTimeStep(Nt, dt);
    
    % compare with Green's function solution
    D1  = thermal_conductivity_1 / (density_1 * specific_heat_1);
    P1  = blood_density * blood_specific_heat * blood_perfusion_rate_1 / (density_1 * specific_heat_1);
    T_exact1 = bioheatExact(37 + s1, 0, [D1, P1, blood_ambient_temp], kgrid.dx, Nt*dt);
    
    D2  = thermal_conductivity_2 / (density_2 * specific_heat_2);
    P2  = blood_density * blood_specific_heat * blood_perfusion_rate_2 / (density_2 * specific_heat_2);
    T_exact2 = bioheatExact(37 + s2, 0, [D2, P2, blood_ambient_temp], kgrid.dx, Nt*dt);
    
    T_exact = zeros(Nx, Ny, Nz);
    T_exact(:, 1:Ny/2, :)     = T_exact1(:, 1:Ny/2, :);
    T_exact(:, Ny/2+1:end, :) = T_exact2(:, Ny/2+1:end, :);
    
    % calculate error
    error_max = abs(kdiff_max.T - T_exact);
    error_min = abs(kdiff_min.T - T_exact);
    
    % error using the exact solution for each side
    error_best = [error_max(:, 1:Ny/2, :), error_min(:, Ny/2+1:end, :)];
    
    % compute the maximum error
    L_inf = max(error_best(:));
    
    % compute pass
    if (L_inf > COMPARISON_THRESH)
        test_pass = false;
    end
    
    % =====================================================================
    % PLOT COMPARISONS
    % =====================================================================
    
    if plot_comparisons
    
        figure;
        
        subplot(2, 3, 1);
        imagesc(squeeze(kdiff_max.T(:, :, Nz/2)));
        axis image;
        colorbar;
        title('k-Wave with Dref max (matches left side)');
    
        subplot(2, 3, 2);
        imagesc(squeeze(kdiff_min.T(:, :, Nz/2)));
        axis image;
        colorbar;
        title('k-Wave with Dref min (matches right side)');    
        
        subplot(2, 3, 3);
        imagesc(squeeze(T_exact(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Exact solution');
    
        subplot(2, 3, 4);
        imagesc(squeeze(error_max(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with Dref max');
    
        subplot(2, 3, 5);
        imagesc(squeeze(error_min(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with Dref min');
        
        subplot(2, 3, 6);
        imagesc(squeeze(error_best(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with correct ref for each half');  
        
    end
    
    % =====================================================================
    % SIMULATION - SOURCE
    % =====================================================================
    
    % set Gaussian heat source
    x_offset = 15e-3;
    y_offset = 15e-3;
    z_offset = 0;
    width = 3*dx;
    source.Q = 3e4 * exp(-( ((kgrid.x-x_offset)/width).^2 + ((kgrid.y-y_offset)/width).^2 + ((kgrid.z-z_offset)/width).^2 ));
    
    % create kWaveDiffusion object and take time steps
    medium.diffusion_coeff_ref = 'max';
    medium.perfusion_coeff_ref = 'min';
    kdiff_max = kWaveDiffusion(kgrid, medium, source, [], input_args{:});
    kdiff_max.takeTimeStep(Nt, dt);
    
    % create kWaveDiffusion object again using opposite reference value and
    % take time steps 
    medium.diffusion_coeff_ref = 'min';
    medium.perfusion_coeff_ref = 'max';
    kdiff_min = kWaveDiffusion(kgrid, medium, source, [], input_args{:});
    kdiff_min.takeTimeStep(Nt, dt);
    
    % compare with Green's function solution
    S = source.Q ./ (medium.density .* medium.specific_heat);
    T_exact1 = bioheatExact(37 + s1, S, [D1, P1, blood_ambient_temp], kgrid.dx, Nt*dt);
    T_exact2 = bioheatExact(37 + s2, S, [D2, P2, blood_ambient_temp], kgrid.dx, Nt*dt);
    
    T_exact = zeros(Nx, Ny, Nz);
    T_exact(:, 1:Ny/2, :)     = T_exact1(:, 1:Ny/2, :);
    T_exact(:, Ny/2+1:end, :) = T_exact2(:, Ny/2+1:end, :);
    
    % calculate error
    error_max = abs(kdiff_max.T - T_exact);
    error_min = abs(kdiff_min.T - T_exact);
    
    % error using the exact solution for each side
    error_best = [error_max(:, 1:Ny/2, :), error_min(:, Ny/2+1:end, :)];
    
    % compute the maximum error
    L_inf = max(error_best(:));
    
    % compute pass
    if (L_inf > COMPARISON_THRESH)
        test_pass = false;
    end
    
    % =====================================================================
    % PLOT COMPARISONS
    % =====================================================================
    
    if plot_comparisons
    
        figure;
        
        subplot(2, 3, 1);
        imagesc(squeeze(kdiff_max.T(:, :, Nz/2)));
        axis image;
        colorbar;
        title('k-Wave with Dref max (matches left side)');
    
        subplot(2, 3, 2);
        imagesc(squeeze(kdiff_min.T(:, :, Nz/2)));
        axis image;
        colorbar;
        title('k-Wave with Dref min (matches right side)');    
        
        subplot(2, 3, 3);
        imagesc(squeeze(T_exact(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Exact solution');
    
        subplot(2, 3, 4);
        imagesc(squeeze(error_max(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with Dref max');
    
        subplot(2, 3, 5);
        imagesc(squeeze(error_min(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with Dref min');
        
        subplot(2, 3, 6);
        imagesc(squeeze(error_best(:, :, Nz/2)));
        axis image;
        colorbar;
        title('Error with correct ref for each half');  
        
    end

end