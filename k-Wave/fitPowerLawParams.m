function [a0_fit, y_fit] = fitPowerLawParams(a0, y, c0, f_min, f_max, plot_fit)
%FITPOWERLAWPARAMS Fit power law absorption parameters for highly absorbing media.
%
% DESCRIPTION:
%     fitPowerLawParams calculates the absorption parameters that should
%     be defined in the simulation functions given the desired power law
%     absorption behaviour defined by a0 and y. This takes into account the
%     actual absorption behaviour exhibited by the fractional Laplacian
%     wave equation. 
%
%     This fitting is required when using large absorption values or high
%     frequencies, as the fractional Laplacian wave equation solved in
%     kspaceFirstOrderND and kspaceSecondOrder no longer encapsulates
%     absorption of the form a = a0*f^y.
%
%     The returned values should be used to define the medium.alpha_coeff
%     and medium.alpha_power within the simulation functions. The
%     absorption behaviour over the frequency range f_min:f_max will then
%     follow the power law defined by a0 and y.
%
% USAGE:
%     [a0_fit, y_fit] = fitPowerLawParams(a0, y, c0, f_min, f_max)
%     [a0_fit, y_fit] = fitPowerLawParams(a0, y, c0, f_min, f_max, plot_fit)
%
% INPUTS:
%     a0          - desired power law absorption prefactor [dB/(MHz^y cm)]
%     y           - desired power law exponent
%     c0          - medium sound speed [m/s]
%     f_min       - minimum frequency of interest [Hz]
%     f_max       - maximum frequency of interest [Hz]
%
% OPTIONAL INPUTS
%     plot_fit    - Boolean controlling whether the final fit is
%                   displayed (default = false)
%
% OUTPUTS:
%     a0_fit      - power law absorption prefactor that should be used to
%                   define medium.alpha_coeff in the simulation functions
%     y_fit       - power law exponent that should be used to define
%                   medium.alpha_power in the simulation functions
%
% ABOUT:
%     author      - Bradley E. Treeby
%     date        - 23rd June 2014
%     last update - 13th January 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2019 Bradley Treeby

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

% check for plot_fit is defined
if nargin < 6
    plot_fit = false;
end

% define frequency axis
f = getSpacedPoints(f_min, f_max, 200);
w = 2 * pi * f;

% convert user defined a0 to Nepers/((rad/s)^y m)
a0_np = db2neper(a0, y);

% define desired absorption behaviour
desired_absorption = a0_np .* w .^ y;

% find the corresponding values of a0 and y that should be used in the
% fractional Laplacian wave equation to give the desired absorption
% behaviour taking into account second order effects
vals = fminsearch(@abs_func, [a0_np, y]);

% assign outputs
a0_np_fit = vals(1);
y_fit = vals(2);

% convert absorption prefactor back to dB/(MHz^y cm)
a0_fit = neper2db(a0_np_fit, y_fit);

% plot the final fit if desired
if plot_fit
    
    % compute absorption behaviour
    absorption_wo_fit = a0_np     * w.^y     ./ ( 1 - (y     + 1) .* a0_np     * c0 * tan(pi * y     / 2) * w.^(y     - 1) );
    absorption_fit    = a0_np_fit * w.^y_fit ./ ( 1 - (y_fit + 1) .* a0_np_fit * c0 * tan(pi * y_fit / 2) * w.^(y_fit - 1) );
    
    % convert from Np/m to dB/cm
    conv_factor = (0.01 * 20 * log10(exp(1)));
    desired_absorption = desired_absorption .* conv_factor;
    absorption_wo_fit  = absorption_wo_fit  .* conv_factor;
    absorption_fit     = absorption_fit     .* conv_factor;
    
    % get suitable x-axis scale factor
    [~, scale, prefix] = scaleSI(f(end));
    
    % downsample factor for plotting actual absorption
    ds = 5;
    
    % plot
    figure;
    plot(f .* scale, desired_absorption, 'k-');
    hold on;
    plot(f .* scale, absorption_wo_fit, 'b--');
    plot(f(1:ds:end).*scale, absorption_fit(1:ds:end), 'rx');
    set(gca, 'FontSize', 12);
    xlabel(['Frequency [' prefix 'Hz]']);
    ylabel('Absorption [dB/cm]');
    legend('Desired Power Law Absorption', 'Actual Absorption Using Original Values', 'Actual Absorption Using Fitted Values', 'Location', 'NorthWest');
    title(['Fitted Absorption Values: \alpha_0 = ' num2str(a0_fit) ', y = ' num2str(y_fit)]);
    axis tight;
    
end


% nested function containing second-order absorption behaviour
% ------------------------------------------------------------
function absorption_error = abs_func(trial_vals)

    % give names to input values
    a0_np_trial = trial_vals(1);
    y_trial  = trial_vals(2);
    
    % evaluate absorption curve based on Eq. (40) in fractional KV paper
    actual_absorption = a0_np_trial * w .^ y_trial ./ ( 1 - (y_trial + 1) .* a0_np_trial * c0 * tan(pi * y_trial / 2) * w.^(y_trial - 1) );
    
    % compute L2 error norm
    absorption_error = sqrt(sum((desired_absorption - actual_absorption).^2));
    
end
% ------------------------------------------------------------

% end of function
end
