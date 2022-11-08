function p_axial = focusedAnnulusONeil(radius, diameter, amplitude, phase, frequency, sound_speed, density, axial_position)
%FOCUSEDANNULUSONEIL Compute axial pressure for focused annulus transducer using O'Neil's solution
%
% DESCRIPTION:
%     focusedAnnulusONeil calculates the axial pressure for a focused
%     annulus transducer using O'Neil's solution (O'Neil, H. Theory of
%     focusing radiators. J. Acoust. Soc. Am., 21(5), 516-526, 1949). The
%     annuluar elements are uniformly driven by a continuous wave sinusoid
%     at a given frequency and normal surface velocity.
%
%     The solution is evaluated at the positions along the beam axis given
%     by axial_position (where 0 corresponds to the transducer surface).
%
%     Note, O'Neil's formulae are derived under the assumptions of the
%     Rayleigh integral, which are valid when the transducer diameter is
%     large compared to both the transducer height and the acoustic
%     wavelength.
%
% USAGE:
%     p_axial = focusedAnnulusONeil(radius, diameter, amplitude, phase, frequency, sound_speed, density, axial_position)
%
% INPUTS:
%     radius           - transducer radius of curvature [m]
%     diameter         - 2 x num_elements array containing pairs of inner
%                        and outer aperture diameter (diameter of opening)
%                        [m].
%     amplitude        - array containing the normal surface velocities for
%                        each element [m/s] 
%     phase            - array containing the phase for each element [rad] 
%     frequency        - driving frequency [Hz]
%     sound_speed      - speed of sound in the propagating medium [m/s]
%     density          - density in the propagating medium [kg/m^3]
%     axial_position   - vector of positions along the beam axis where the
%                        pressure amplitude is calculated [m]
%
% OUTPUTS:
%     p_axial          - pressure amplitude at the positions specified by
%                        axial_position [Pa]
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 1st April 2021
%     last update      - 3rd August 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2021 Bradley Treeby
%
% See also focusedBowlONeil.

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

% make sure diameter input is the correct dimension
validateattributes(diameter, {'numeric'}, {'nonnegative', 'real', 'finite', 'size', [2, NaN]}, mfilename, 'diameter');

% make sure amplitude and phase have the correct length for the number of
% elements
num_elements = size(diameter, 2);
validateattributes(amplitude, {'numeric'}, {'nonnegative', 'real', 'finite', 'numel', num_elements}, mfilename, 'amplitude');
validateattributes(phase, {'numeric'}, {'nonnegative', 'real', 'finite', 'numel', num_elements}, mfilename, 'phase');

% pre-allocate output
p_axial = zeros(size(axial_position));

% loop over elements and sum fields
for ind = 1:num_elements
    
    % get complex pressure for bowls with inner and outer aperature diam
    if diameter(1, ind) == 0
        p_el_inner = 0;
    else
        [~, ~, p_el_inner] = focusedBowlONeil(radius, diameter(1, ind), amplitude(ind), frequency, sound_speed, density, axial_position, []);
    end
    [~, ~, p_el_outer] = focusedBowlONeil(radius, diameter(2, ind), amplitude(ind), frequency, sound_speed, density, axial_position, []);
    
    % pressure for annular element
    p_el = p_el_outer - p_el_inner;
    
    % account for phase
    p_el = abs(p_el) .* exp(1i .* (angle(p_el) + phase(ind)));
    
    % add to complete response
    p_axial = p_axial + p_el;
    
end

% take magnitude
p_axial = abs(p_axial);
