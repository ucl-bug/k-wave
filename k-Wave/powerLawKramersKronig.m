function c_kk = powerLawKramersKronig(w, w0, c0, a0, y)
%POWERLAWKRAMERSKRONIG Calculate dispersion for power law absorption.
%
% DESCRIPTION:
%     powerLawKramersKronig computes the variation in sound speed for an
%     attenuating medium using the Kramers-Kronig for power law
%     attenuation where att = a0*w^y. The power law parameters must be in
%     Nepers / m, with the frequency in rad/s. The variation is given
%     about the sound speed c0 at a reference frequency w0.
%
% USAGE:
%     c_kk = powerLawKramersKronig(w, w0, c0, a0, y)
%
% INPUTS:
%     w             - input frequency array [rad/s]
%     w0            - reference frequency [rad/s]
%     c0            - sound speed at w0 [m/s]
%     a0            - power law coefficient [Nepers/((rad/s)^y m)]
%     y             - power law exponent, where 0 < y < 3
%
% OUTPUTS:
%     c_kk          - variation of sound speed with w [m/s]
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 8th August 2008
%     last update   - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2008-2017 Bradley Treeby
%
% See also db2neper

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

if y >= 3 || y <= 0
    
    % for y values outside range, assign c0 to output
    disp('WARNING: y must be 0 < y < 3.');
    c_kk = c0 * ones(size(w));
    
elseif y == 1
    
    % Kramers-Kronig for y = 1
    c_kk = 1 ./ ( 1 ./ c0 - 2 .* a0 .* log(w ./ w0) ./ pi );
    
else
    
    % Kramers-Kronig for 0 < y < 1 and 1 < y < 3
    c_kk = 1 ./( 1 ./ c0 + a0 .* tan(y .* pi ./ 2) .* (w.^(y - 1) - w0.^(y - 1)) );
    
end