function p_mendousse = mendousse(x, t, source_freq, p0, c0, rho0, BonA, alpha_0)
%MENDOUSSE Compute Mendousse's solution for nonlinear wave propagation in viscous media.  
%
% DESCRIPTION:
%     mendousse calculates the propagation of a monofrequency plane wave
%     source in a thermoviscous medium with absorption given by
%     alpha_0*f^2. The solution is taken from Eq. (264) in Chapter 4 of
%     Nonlinear acoustics by Hamilton and Blackstock (2008). The infinite
%     sum is adaptively truncated when the moving average of the previous
%     five sum terms is smaller than a predefined convergence percentage
%     (0.01 percent by default).
%
% USAGE:
%     p_mendousse = mendousse(x, t, source_freq, p0, c0, rho0, BonA, alpha_0)
%
% INPUTS:
%     x           - position [m]
%     t           - time [s]
%     source_freq - frequency of plane wave [Hz]
%     p0          - source pressure [Pa]
%     c0          - medium sound speed [m/s]
%     rho0        - medium density [kg/m^3]
%     BonA        - nonlinearity parameter B/A
%     alpha_0     - absorption coefficient [dB/(MHz^2 cm)]
%
% OUTPUTS:
%     p_mendousse - calculated pressure field
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 21st February 2011
%     last update - 7th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby
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

% stop the summation when the current term contributes less than this amount
CONVERGENCE_PERCENTAGE = 0.01;   

% minimum number of sum terms to use
MIN_NUM_SUM_TERMS = 20;

% number of pressure values to use in the moving average
MOVING_AV = 5;

% calculate the mach number
mach_num = p0 ./ (rho0 .* c0^2);

% calculate the nonlinearity coefficient
beta = 1 + BonA ./ 2;

% convert the absorption to nepers per meter
alpha_0 = db2neper(alpha_0, 2);

% calculate the absorption parameter
alpha = alpha_0 .* (2 .* pi .* source_freq).^2;

% calculate the Goldberg number
goldberg_num = beta .* mach_num .* 2 .* pi .* source_freq ./ (c0 .* alpha);

% preallocate output variable
p_mendousse = zeros(1, length(x));

% calculate series
for loop_index = 1:length(x)
   
    % =====================================================================
    % COMPUTE NUMERATOR SUMMATION
    % =====================================================================
    
    % initialise loop variables
    p_term1 = 0;                    % pressure
    sum_term_contrib = 1;           % percentage contribution of current sum term
    n = 1;                          % summation index variable
    p_prev = zeros(1, MOVING_AV);   % moving average variable

    % loop the sum until it reaches an acceptable level of convergence
    while (sum_term_contrib > CONVERGENCE_PERCENTAGE) || (n < MIN_NUM_SUM_TERMS)

        % compute next sum term and add to the total
        p_term1 = p_term1 + (-1)^(n + 1) .* n .* besseli(n, goldberg_num ./ 2) ... 
            .* exp(-n.^2 .* alpha .* x(loop_index)) ...
            .* sin(n .* 2 .* pi .* source_freq .* t(loop_index));
        
        % update moving average of the contribution of sum term in a fifo
        % paradigm
        p_prev = circshift(p_prev, [0, 1]);
        p_prev(1) = p_term1;
        sum_term_contrib = 100 .* abs((p_term1 - mean(p_prev)) ./ p_term1);
        
        % update summation variable
        n = n + 1;

    end
    
    % =====================================================================
    % COMPUTE DENOMINATOR SUMMATION
    % =====================================================================
    
    % initialise loop variables
    p_term2 = 0;                    % pressure
    sum_term_contrib = 1;           % percentage contribution of current sum term
    n = 1;                          % summation index variable
    p_prev = zeros(1, MOVING_AV);   % moving average variable

    % loop the sum until it reaches an acceptable level of convergence
    while (sum_term_contrib > CONVERGENCE_PERCENTAGE) || (n < MIN_NUM_SUM_TERMS)

        % compute next sum term and add to the total
        p_term2 = p_term2 + (-1)^(n) * besseli(n, goldberg_num ./ 2) ...
            .* exp(-n^2 .* alpha .* x(loop_index)) ...
            .* cos(n .* 2 .* pi .* source_freq .* t(loop_index)); 
        
        % update moving average of the contribution of sum term in a fifo
        % paradigm
        p_prev = circshift(p_prev, [0, 1]);
        p_prev(1) = p_term2;
        sum_term_contrib = 100 .* abs((p_term2 - mean(p_prev)) ./ p_term2);
        
        % update summation variable
        n = n + 1;

    end
     
    % =====================================================================
    % APPLY SCALING PARAMETERS AND STORE THE VALUE OF P
    % =====================================================================
    
    % check if the numerator is zero
    if p_term1 ~= 0
        p_term1 = 4 .* p_term1 ./ goldberg_num;
    else
        p_term1 = 0;
    end
    p_mendousse(loop_index) = p_term1 ./ (besseli(0, goldberg_num ./ 2) + 2 .* p_term2);
    
end

% scale output by p0
p_mendousse = p0 .* p_mendousse;