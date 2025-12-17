function test_pass = angularSpectrumCW_absorption(~, ~)
% DESCRIPTION:
%     Unit test to propagate a plane wave through an absorbing medium and
%     check that the absorption of the wave amplitude agrees with the
%     specified absorption value.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 6th March 2018
%     last update - 19th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018- Bradley Treeby

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

% set comparison threshold
COMPARISON_THRESH = 1e-14;

% set pass variable
test_pass = true;

% simulation settings
Nx = 64;
Ny = 64;
Nz = 12;
dx = 1e-4;
f0 = 2e6;

% assign medium properties
medium.sound_speed = 1500;
medium.alpha_coeff = 0.75;
medium.alpha_power = 1.5;

% generate input plane
input_plane = ones(Nx, Ny);

% run simulation
pressure = angularSpectrumCW(input_plane, dx, (0:(Nz - 1)) * dx, f0, medium, 'FFTLength', Nx);

% get output pressure magnitude
p = abs(pressure(Nx/2, Ny/2, Nz));

% desired absorption in dB/cm
a_des = medium.alpha_coeff .* (f0 * 1e-6) .^ medium.alpha_power;

% actual absorption in dB/cm
dist = (Nz - 1) * dx * 1e2;
a_act = -20 * log10(p) / dist;

% compute error
err = abs(a_des - a_act) / abs(a_des);

% check for test pass
if err > COMPARISON_THRESH
    test_pass = false;
end