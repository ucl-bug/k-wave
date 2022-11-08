function BonA = waterNonlinearity(T)
%WATERNONLINEARITY Calculate B/A of water with temperature.
%
% DESCRIPTION:
%     waterNonlinearity calculates the parameter of nonlinearity B/A at a
%     given temperature using a fourth-order polynomial fitted to the data
%     given by Beyer (1960).
%
% USAGE:
%     BonA = waterNonlinearity(T)
%
% INPUTS:
%     T             - water temperature in the range 0 to 100 [degC]
%
% OUTPUTS:
%     BonA          - parameter of nonlinearity
%
% ABOUT:
%     author        - Bradley E. Treeby
%     date          - 22nd February 2018
%     last update   - 4th April 2019
%
% REFERENCES:
%     [1] R. T Beyer (1960) "Parameter of nonlinearity in fluids," J.
%     Acoust. Soc. Am., 32(6), 719-721.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby
%
% See also waterAbsorption, waterDensity, waterSoundSpeed

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
if any(T(:) > 100) || any(T(:) < 0)
    error('T must be between 0 and 100.');
end

% find value
p(1) = -4.587913769504693e-08;
p(2) = 1.047843302423604e-05;
p(3) = -9.355518377254833e-04;
p(4) = 5.380874771364909e-2;
p(5) = 4.186533937275504;
BonA = polyval(p, T);

% % original data
% T = [0, 10, 20, 30, 40, 60, 80, 100];
% BonA = [4.2, 4.6, 5.0, 5.2, 5.4, 5.7, 6.0, 6.1];
% 
% % fitting (requires curve fitting toolbox)
% f0 = fit(T.', BonA.', 'poly4');
% 
% % plot fit and original data
% figure;
% plot(T, BonA, 'rx');
% hold on;
% T_interp = 0:0.1:100;
% BonA_interp = polyval([f0.p1, f0.p2, f0.p3, f0.p4, f0.p5], T_interp);
% plot(T_interp, BonA_interp, 'k-');
% xlabel('Temperature [degC]');
% ylabel('BonA');