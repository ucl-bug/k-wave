function density = waterDensity(T)
%WATERDENSITY Calculate density of air-saturated water with temperature.
%
% DESCRIPTION:
%     waterDensity calculates the density of air-saturated water at a given
%     temperature using the 4th order polynomial given by Jones (1992).
%
% USAGE:
%     density = waterDensity(T)
%
% INPUTS:
%     T             - water temperature in the range 5 to 40 [degC]
%
% OUTPUTS:
%     density       - density of water [kg/m^3]
%
% ABOUT:
%     author        - Bradley E. Treeby
%     date          - 22nd February 2018
%     last update   - 4th April 2019
%
% REFERENCES:
%     [1] F. E. Jones and G. L. Harris (1992) "ITS-90 Density of Water
%     Formulation for Volumetric Standards Calibration," J. Res. Natl.
%     Inst. Stand. Technol., 97(3), 335-340.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby
%
% See also waterAttenuation, waterNonlinearity, waterSoundSpeed

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

% check limits
if any(T(:) > 40) || any(T(:) < 5)
    error('T must be between 5 and 40.');
end

% calculate density of air-saturated water
density = 999.84847 + 6.337563e-2 .* T ...
                    - 8.523829e-3 .* T.^2 ...
                    + 6.943248e-5 .* T.^3 ...
                    - 3.821216e-7 .* T.^4;