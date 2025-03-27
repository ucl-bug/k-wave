function test_pass = kWaveDiffusion_compare_with_pde(~, plot_simulations)
% DESCRIPTION:
%       Unit test to compare kWaveDiffusion with a solution to the heat
%       equation solved using the MATLAB PDE toolbox.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 24 July 2020
%       last update - 24 July 2020
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2020-2020 Bradley Treeby

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

% check for PDE toolbox
v = ver;
if ~any(strcmp('Partial Differential Equation Toolbox', {v.Name}))
    warning('kWaveDiffusion_compare_with_pde not tested. The MATLAB Partial Differential Equation Toolbox must be installed to run this test.');
    return
end

% check for dtt functions
if ~exist('dtt1D', 'file') 
    warning('kWaveDiffusion_compare_with_pde not tested. The MATLAB dtt library must be installed to run this test.');
    return
end

% set comparison threshold
COMPARISON_THRESH = 0.01;

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_simulations = true;
end

% =========================================================================
% DEFINE PROBLEM
% =========================================================================

% domain size
x_size = 0.1;
y_size = 0.1;

% mesh size (dx for k-Wave, maximum element size for MATAB)
dx = 5e-3;

% time
Nt = 200;
dt = 5;

% material properties
thermal_conductivity = 0.52;
density              = 1079;
specific_heat        = 3540;

% specify initial temperature as a Gaussian (use anonymous function)
width = 4 * dx;
T0 = @(location) exp( -(location.x ./ width).^2 - (location.y ./width).^2 );

% =========================================================================
% kWAVEDIFFUSION SOLUTION
% =========================================================================

% compute grid size
Nx = round(x_size / dx);
Ny = round(y_size / dx);

% create grid structure
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% define time array
t_array = (0:(Nt - 1)) * dt;

% create initial temperature distribution
source.T0 = T0(kgrid);

% plot the initial temperature
if plot_simulations
    figure;
    imagesc(kgrid.y_vec, kgrid.x_vec, source.T0);
    axis image;
    colorbar;
    title('Initial temperature');
end

% define the medium properties
medium.thermal_conductivity = thermal_conductivity;
medium.density              = density;
medium.specific_heat        = specific_heat;

% define boundary condition
medium.boundary_condition = 'conducting';

% define the sensor
sensor.mask = ones(Nx, Ny);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, sensor, ...
    'PlotSim', false);

% take time steps
kdiff.takeTimeStep(Nt, dt);

% reshape the temperature
T_kwave = reshape(kdiff.sensor_data, Nx, Ny, Nt);

% =========================================================================
% FEM SOLUTION
% =========================================================================

% create model
thermalmodel = createpde('thermal', 'transient');

% define geometry using constructive solid geometry: rectangle is code 3, 4
% sides, followed by x-coordinates and then y-coordinates 
% make the FEM model one grid point bigger, as  WAWA symmetry has the BC
% outside
geom = [3, 4, ...
    kgrid.x_vec(1) - dx, kgrid.x_vec(end) + dx, kgrid.x_vec(end) + dx, kgrid.x_vec(1) - dx, ...
    kgrid.y_vec(1) - dx, kgrid.y_vec(1) - dx,   kgrid.y_vec(end) + dx, kgrid.y_vec(end) + dx]';
geom = decsg(geom);

% create rectangular geometry from edges
geometryFromEdges(thermalmodel, geom);

% plot the geometry
if plot_simulations
    figure;
    pdegplot(thermalmodel, 'EdgeLabels', 'on');
    axis image;
end

% specifty thermal properties
thermalProperties(thermalmodel, ...
    'ThermalConductivity', thermal_conductivity, ...
    'MassDensity', density, ...
    'SpecificHeat', specific_heat);

% specify internal heat sourt                           
internalHeatSource(thermalmodel, 0);

% specify initial temperature
thermalIC(thermalmodel, T0);

% define Dirichlet boundary conditions
% Neumann -> gradient zero -> insulating -> radiation BC?
% Dirichlet -> value is zero -> conducting -> temperature BC?
thermalBC(thermalmodel, 'Edge', 1:4, 'Temperature', 0);

% generate mesh with no element larger than the dx used in kWaveDiffusion
generateMesh(thermalmodel, 'Hmax', dx);

% plot the mesh
if plot_simulations
    figure;
    pdemesh(thermalmodel); 
    axis equal;
end

% generate the transient solution
results = solve(thermalmodel, t_array);

% % plot the field
% T = results.Temperature;
% figure;
% for ind = 1:10:Nt
%     pdeplot(thermalmodel, 'XYData', T(:, ind));
%     axis image;
%     title(['Time = ' scaleTime(t_array(ind))]);
%     drawnow;
%     pause(0.2);
% end

% interpolate the solution onto a regular mesh
T_fem = interpolateTemperature(results, kgrid.x, kgrid.y, 1:Nt);
T_fem = reshape(T_fem, Nx, Ny, Nt);

% compute error
err = max(abs(T_fem(:) - T_kwave(:)));

if err > COMPARISON_THRESH
    test_pass = false;
end

% =========================================================================
% PLOTTING
% =========================================================================

if plot_simulations
    
    % plot the field
    figure;
    for ind = 1:10:Nt
        
        subplot(2, 2, 1);
        imagesc(kgrid.y_vec, kgrid.x_vec, T_fem(:, :, ind));
        axis image;
        title('FEM');
        colorbar;
        
        subplot(2, 2, 2);
        imagesc(kgrid.y_vec, kgrid.x_vec, T_kwave(:, :, ind));
        axis image;
        title('k-Wave');
        colorbar;    
        
        subplot(2, 2, 3);
        imagesc(kgrid.y_vec, kgrid.x_vec, abs(T_fem(:, :, ind) - T_kwave(:, :, ind)));
        axis image;
        title('Difference');
        colorbar;
        
        subplot(2, 2, 4);
        plot(kgrid.x_vec, squeeze(T_fem(:, Ny/2, ind)), 'k-');
        hold on;
        plot(kgrid.x_vec, squeeze(T_kwave(:, Ny/2, ind)), 'bo', 'MarkerSize', 12);
        legend('FEM', 'k-Wave');
        hold off;
        set(gca, 'YLim', [0, 1]);
        
        drawnow;
        pause(0.2);
        
    end

end