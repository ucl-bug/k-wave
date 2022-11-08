function T = bioheatExact(T0, S, material, dx, t)
%BIOHEATEXACT Compute exact solution to Pennes' bioheat equation in homogeneous media. 
%
% DESCRIPTION:
%     bioheatExact calculates the exact solution to Pennes' bioheat
%     equation in a homogeneous medium on a uniform Cartesian grid using a
%     Fourier-based Green's function solution assuming a periodic boundary
%     condition [1]. The function supports inputs in 1D, 2D, and 3D. The
%     exact equation solved is given by  
%
%         dT/dt = D * d^2T/dx^2 - P * (T - Ta) + S
%
%     where the coefficients are defined below. Pennes' bioheat equation is
%     often given in the alternative form 
%
%         P0 * C0 * dT/dt =  Kt * d^2T/dx^2 - Pb * Wb * Cb * (T - Ta) + Q
%
%         T:  temperature                     [degC]
%         C0: tissue specific heat capacity   [J/(kg.K)]
%         P0: tissue density                  [kg/m^3] 
%         Kt: tissue thermal conductivity     [W/(m.K)]
%         Pb: blood density                   [kg/m^3] 
%         Wb: blood perfusion rate            [1/s]
%         Ta: blood arterial temperature      [degC]
%         Cb: blood specific heat capacity    [J/(kg.K)]
%         Q:  volume rate of heat deposition  [W/m^3]
%
%     In this case, the function inputs are calculated by
%
%         D = Kt / (P0 * C0);
%         P = Pb * Wb * Cb / (P0 * C0);
%         S = Q / (P0 * C0);
%
%     If the perfusion coefficient P is set to zero, bioheatExact
%     calculates the exact solution to the heat equation in a homogeneous
%     medium.
%
%     [1] Gao, B., Langer, S., & Corry, P. M. (1995). Application of the
%     time-dependent Green's function and Fourier transforms to the
%     solution of the bioheat equation. International Journal of
%     Hyperthermia, 11(2), 267-285.
%
% USAGE:
%     T = bioheatExact(T0, S, material, dx, t)
%     T = bioheatExact(T0, S, [D, P, Ta], dx, t)
%
% INPUTS:
%     T0          - matrix of the initial temperature distribution at each
%                   grid point [degC] 
%     S           - matrix of the heat source at each grid point [K/s]
%     material    - material coefficients given as a three element vector
%                   in the form: material = [D, P, Ta], where
%
%                       D:  diffusion coefficient [m^2/s]
%                       P:  perfusion coefficient [1/s]
%                       Ta: arterial temperature  [degC]
%
%     dx          - grid point spacing [m]
%     t           - time at which to calculate the temperature field [s]
%
% OUTPUTS:
%     T           - temperature field at time t [degC]
%
% ABOUT:
%     author      - Bradley Treeby and Teedah Saratoon
%     date        - 18th August 2015
%     last update - 11th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Bradley Treeby and Teedah Saratoon
%
% See also kWaveDiffusion

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

% check that T0 is a matrix
if numel(T0) == 1
    error('T0 must be defined as a matrix.');
end

% if S is not 0, check that S and T0 are the same size
if ~((numel(S) == 1) && (S == 0)) && ~all(size(S) == size(T0))
    error('T0 and S must be the same size.')
end

% extract material properties
D  = material(1);   % diffusion coefficient
P  = material(2);   % perfusion coefficient
Ta = material(3);   % blood arterial temperature

% check the medium properties are homogeneous
if numel(P) > 1 || numel(D) > 1 || numel(Ta) > 1
    error('Medium properties must be homogeneous.');
end

% create the set of wavenumbers
kgrid = kWaveGrid(size(T0, 1), dx, size(T0, 2), dx, size(T0, 3), dx);

% define Green's function propagators
T0_propagator = exp(-(D .* ifftshift(kgrid.k).^2 + P) .* t);
Q_propagator  = (1 - T0_propagator) ./ (D .* ifftshift(kgrid.k).^2 + P);

% replace Q propagator with limits for k == 0
%     if P == 0, the limit is t
%     if P ~= 0, the limit is (1 - exp(-P*t))/P
if (numel(P) == 1) && (P == 0)
    Q_propagator(isnan(Q_propagator)) = t;
else
    Q_propagator(isnan(Q_propagator)) = (1 - exp(-P .* t)) ./ P;
end

% calculate exact Green's function solution (Eq. 12 [1])
if (numel(S) == 1) && (S == 0)
    T_change = real(ifftn( T0_propagator .* fftn(T0 - Ta) ));    
else
    T_change = real(ifftn( T0_propagator .* fftn(T0 - Ta) + Q_propagator .* fftn(S) ));
end
T = T_change + Ta;