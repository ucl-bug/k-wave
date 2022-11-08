% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to modify
%     p_plot if optional input 'PlotScale' is set to 'auto' or 'decibel'.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 4th February 2011
%     last update - 15th May 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2018 Bradley Treeby

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

% update plot scale if set to automatic, setting the plot scale to be
% symmetric about 0
if flags.plot_scale_auto
    if ~flags.elastic_code
        
        % get maximum value
        mx = max(abs(p_plot(:)));
        
        % offset scale if maximum is zero
        if mx == 0
            mx = 0.1;
        end
        
        % set plot scale
        plot_scale = [-mx, mx];
        
    else
        
        % get maximum value
        mx_ii = max(abs(sii_plot(:)));
        mx_ij = max(abs(sij_plot(:)));
        
        % offset scale if maximum is zero
        if mx_ii == 0
            mx_ii = 0.1;
        end
        if mx_ij == 0
            mx_ij = 0.1;
        end
        
        % set plot scale
        plot_scale = [-mx_ii, mx_ii, -mx_ij, mx_ij];
        
    end
end

% rescale plot variable if plot scale set to log
if flags.plot_scale_log
    
    % update plot scale if set to auto
    if flags.plot_scale_auto
        alt_plot_scale_lin = plot_scale;
        alt_plot_scale_log = log10(abs(plot_scale) + log_scale_comp_factor) - log10(log_scale_comp_factor);
        alt_plot_scale_log(1) = -alt_plot_scale_log(1);
    else
        
        % truncate data to the given plot scale
        p_plot(p_plot < alt_plot_scale_lin(1)) = alt_plot_scale_lin(1);
        p_plot(p_plot > alt_plot_scale_lin(2)) = alt_plot_scale_lin(2);
        
    end

    % scale both positive and negative data separately
    p_plot(p_plot > 0) =   log10( p_plot(p_plot > 0) + log_scale_comp_factor) - log10(log_scale_comp_factor);
    p_plot(p_plot < 0) = -(log10(-p_plot(p_plot < 0) + log_scale_comp_factor) - log10(log_scale_comp_factor));
    
    % update plot scale
    plot_scale = alt_plot_scale_log;
    
end