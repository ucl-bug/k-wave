function dt_stability_limit = checkStability(kgrid, medium)
% CHECKSTABILITY Compute maximum stable timestep for k-space fluid models.
%
% DESCRIPTION:
%     checkStability calculates the maximum time step for which the k-space
%     propagation models kspaceFirstOrder1D, kspaceFirstOrder2D and
%     kspaceFirstOrder3D are stable. These models are unconditionally
%     stable when the reference sound speed is equal to or greater than the
%     maximum sound speed in the medium and there is no absorption.
%     However, when the reference sound speed is less than the maximum
%     sound speed the model is only stable for sufficiently small time
%     steps. The criterion is more stringent (the time step is smaller) in
%     the absorbing case.
%
%     The time steps given are accurate when the medium properties are
%     homogeneous. For a heterogeneous media they give a useful, but not
%     exact, estimate.  
%
% USAGE:
%     dt_stability_limit = checkStability(kgrid, medium)
%
% INPUTS:
%     kgrid              - k-Wave grid object return by kWaveGrid
%     medium             - structure containing the medium properties
%
% OUTPUTS:
%     dt_stability_limit - the maximum time step for which the models are
%                          stable. This is set to Inf when the model is
%                          unconditionally stable.
%
% ABOUT:
%     author             - Ben Cox
%     date               - 12th August 2014
%     last update        - 27th February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Ben Cox
%
% See also kspaceFirstOrder1D, kspaceFirstOrder2D, kspaceFirstOrder3D,
% kWaveGrid

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

% define literals
FIXED_POINT_ACCURACY = 1e-12;

% find the maximum wavenumber
kmax = max(kgrid.k(:));

% calculate the reference sound speed for the fluid code, using the
% maximum by default which ensures the model is unconditionally stable
if isfield(medium, 'sound_speed_ref')
    if isnumeric(medium.sound_speed_ref)
        c_ref = medium.sound_speed_ref;
    elseif strcmp(medium.sound_speed_ref, 'min')
        c_ref = min(medium.sound_speed(:));
    elseif strcmp(medium.sound_speed_ref, 'mean')
        c_ref = mean(medium.sound_speed(:));
    elseif strcmp(medium.sound_speed_ref, 'max')
        c_ref = max(medium.sound_speed(:));        
    end
else
    c_ref = max(medium.sound_speed(:));
end

% calculate the timesteps required for stability
if ~isfield(medium, 'alpha_coeff') || all(medium.alpha_coeff(:) == 0)

    % =====================================================================
    % NON-ABSORBING CASE
    % =====================================================================
    
    if c_ref >= max(medium.sound_speed(:))
        
        % set the timestep to Inf when the model is unconditionally stable
        dt_stability_limit = Inf;
        
    else
        
        % set the timestep required for stability when c_ref~=max(medium.sound_speed(:))
        dt_stability_limit = 2/(c_ref * kmax) * asin(c_ref/max(medium.sound_speed(:)));
        
    end
    
else

    % =====================================================================
    % ABSORBING CASE
    % =====================================================================

    % convert the absorption coefficient to nepers.(rad/s)^-y.m^-1
    medium.alpha_coeff = db2neper(medium.alpha_coeff, medium.alpha_power);

    % calculate the absorption constant
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_absorption'))
        absorb_tau = -2 * medium.alpha_coeff .* medium.sound_speed.^(medium.alpha_power - 1);
    else
        absorb_tau = 0;
    end

    % calculate the dispersion constant
    if ~(isfield(medium, 'alpha_mode') && strcmp(medium.alpha_mode, 'no_dispersion'))
        absorb_eta = 2 * medium.alpha_coeff .* medium.sound_speed.^(medium.alpha_power) * tan(pi * medium.alpha_power / 2);
    else
        absorb_eta = 0;
    end
    
    % estimate the timestep required for stability in the absorbing case by
    % assuming the k-space correction factor, kappa = 1 (note that
    % absorb_tau and absorb_eta are negative quantities) 
    temp1 = max(medium.sound_speed(:)) * min(absorb_tau(:)) * kmax^(medium.alpha_power - 1);
    temp2 = 1 - min(absorb_eta(:)) * kmax^(medium.alpha_power - 1);
    dt_estimate = (temp1 + sqrt(temp1^2 + 4 * temp2)) / (temp2 * kmax * max(medium.sound_speed(:)));
    
    % use a fixed point iteration to find the correct timestep, assuming
    % now that kappa = kappa(dt), using the previous estimate as a starting
    % point
    
    % first define the function to iterate
    kappa = @(dt) sinc(c_ref * kmax * dt / 2);
    temp3 = @(dt) max(medium.sound_speed(:)) * min(absorb_tau(:)) * kappa(dt) * kmax^(medium.alpha_power - 1);
    func_to_solve = @(dt) (temp3(dt) + sqrt((temp3(dt))^2 + 4 * temp2 )) / (temp2 * kmax * kappa(dt) * max(medium.sound_speed(:)));

    % run the fixed point iteration  
    dt_stability_limit = dt_estimate;
    dt_old = 0;
    while abs(dt_stability_limit - dt_old) > FIXED_POINT_ACCURACY
        dt_old = dt_stability_limit;
        dt_stability_limit = func_to_solve(dt_stability_limit);
    end
    
end