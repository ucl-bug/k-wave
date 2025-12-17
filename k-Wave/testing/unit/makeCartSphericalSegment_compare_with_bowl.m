function test_pass = makeCartSphericalSegment_compare_with_bowl(plot_comparisons, ~)
%makeCartSphericalSegment_COMPARE_WITH_BOWL Check multiple annuli with bowl.
%
% DESCRIPTION:
%     makeCartSphericalSegment_compare_with_bowl compares the Cartesian
%     points of multiple annuli (spherical segments) with a single bowl.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 30th March 2021
%     last update      - 8th June 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2021- Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set pass variable
test_pass = true;

% set comparison threshold
COMPARISON_THRESH = 5e-4;

% bowl properties
bowl_pos = [0, 0, 0];
focus_pos = [0, 0, 1];
num_points_vec = [69, 70, 200, 1000];
roc = 65e-3;
ap_diam1 = 30e-3;
ap_diam2 = 45e-3;
ap_diam3 = 60e-3;

% loop over some different numbers of points
for points_ind = 1:length(num_points_vec)
    
    num_points = num_points_vec(points_ind);

    % make regular bowl
    bowl_points = makeCartBowl(bowl_pos, roc, ap_diam3, focus_pos, num_points);

    % =====================================================================
    % FULL BOWL
    % =====================================================================

    % make a regular bowl using the annular function
    bowl_annulus_points = makeCartSphericalSegment(bowl_pos, roc, 0, ap_diam3, focus_pos, num_points);

    % check the points are the same
    err = max(abs(bowl_points(:) - bowl_annulus_points(:)));

    % check error
    if err > COMPARISON_THRESH
        test_pass = false;
    end

    % =====================================================================
    % 2-ELEMENT ANNULUS
    % =====================================================================

    % find position where to split the bowl
    bowl_height1 = roc - sqrt(roc^2 - (ap_diam1/2)^2);
    bowl_height2 = roc - sqrt(roc^2 - (ap_diam2/2)^2);

    % split array
    bowl_points_ann1 = bowl_points(:, bowl_points(3, :) <= bowl_height1);
    bowl_points_ann2 = bowl_points(:, (bowl_points(3, :) <=  bowl_height2) & (bowl_points(3, :) > bowl_height1));
    bowl_points_ann3 = bowl_points(:, bowl_points(3, :) >  bowl_height2);

    % count points
    num_points_ann1 = size(bowl_points_ann1, 2);
    num_points_ann2 = size(bowl_points_ann2, 2);
    num_points_ann3 = size(bowl_points_ann3, 2);

    % make the annuli separately
    annulus_points_ann1 = makeCartSphericalSegment(bowl_pos, roc, 0,        ap_diam1, focus_pos, num_points_ann1);
    annulus_points_ann2 = makeCartSphericalSegment(bowl_pos, roc, ap_diam1, ap_diam2, focus_pos, num_points_ann2, [], num_points_ann1);
    annulus_points_ann3 = makeCartSphericalSegment(bowl_pos, roc, ap_diam2, ap_diam3, focus_pos, num_points_ann3, [], num_points_ann1 + num_points_ann2);

    % create the figure
    if plot_comparisons
        figure;
        plot(bowl_points_ann1(1, :), bowl_points_ann1(2, :), 'b.');
        hold on;
        plot(bowl_points_ann2(1, :), bowl_points_ann2(2, :), 'r.');
        plot(bowl_points_ann3(1, :), bowl_points_ann3(2, :), 'k.');
        plot(annulus_points_ann1(1, :), annulus_points_ann1(2, :), 'bo');
        plot(annulus_points_ann2(1, :), annulus_points_ann2(2, :), 'ro');
        plot(annulus_points_ann3(1, :), annulus_points_ann3(2, :), 'ko');
        viscircles([0, 0], ap_diam1/2, 'Color', 'b');
        viscircles([0, 0], ap_diam2/2, 'Color', 'r');
        viscircles([0, 0], ap_diam3/2, 'Color', 'k');
        xlabel('m');
        ylabel('m');
        zlabel('m');
        axis equal;
        grid on;
        box on;
        legend('Bowl Annulus 1', 'Bowl Annulus 2', 'Bowl Annulus 3', 'Annulus 1', 'Annulus 2', 'Annulus 3');
        view([0, 90]);
    end

    % check the points are the same
    err1 = max(abs(bowl_points_ann1(:) - annulus_points_ann1(:)));
    err2 = max(abs(bowl_points_ann2(:) - annulus_points_ann2(:)));
    err3 = max(abs(bowl_points_ann3(:) - annulus_points_ann3(:)));

    % check error
    if (err1 > COMPARISON_THRESH) || (err2 > COMPARISON_THRESH) || (err3 > COMPARISON_THRESH)
        test_pass = false;
    end

end
        