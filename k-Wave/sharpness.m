function s = sharpness(im, metric)
%SHARPNESS Calculate image sharpness metric.
% 
% DESCRIPTION:
%     sharpness returns a scalar metric related to the sharpness of the 2D
%     or 3D image matrix defined by im. By default, the metric is based on
%     the Brenner gradient which returns the sum of the centered
%     finite-difference at each matrix element in each Cartesian direction.
%     Metrics calculated using the Sobel operator or the normalised
%     variance can also be returned by setting the input paramater metric.
%
%     For further details, see B. E. Treeby, T. K. Varslot, E. Z. Zhang,
%     J. G. Laufer, and P. C. Beard, "Automatic sound speed selection in
%     photoacoustic image reconstruction using an autofocus approach," J.
%     Biomed. Opt., vol. 16, no. 9, p. 090501, 2011.
%
% USAGE:
%     s = sharpness(im)
%     s = sharpness(im, metric)
%
% INPUTS:
%     im          - 2D or 3D image data to evaluate
%
% OPTIONAL INPUTS:
%     metric      - sharpness metric. Supported values are: 
%
%                   'Brenner'           (default)
%                   'Tenenbaum'
%                   'NormVariance'
%
% OUTPUTS:
%     s           - computed sharpness metric
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 15th January 2012
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

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

% check for metric input
if nargin == 1
   metric = 'Brenner'; 
end

switch metric
    
    case 'Brenner'
        
        % compute sharpness metric based on the Brenner gradient
        switch numDim(im)
            case 2
                
                % compute metric
                bren_x = (im(1:end - 2, :) - im(3:end, :)).^2;
                bren_y = (im(:, 1:end - 2) - im(:, 3:end)).^2;
                s = sum(bren_x(:)) + sum(bren_y(:));                  
                
            case 3
                
                % compute metric
                bren_x = (im(1:end - 2, :, :) - im(3:end, :, :)).^2;
                bren_y = (im(:, 1:end - 2, :) - im(:, 3:end, :)).^2;
                bren_z = (im(:, :, 1:end - 2) - im(:, :, 3:end)).^2;
                s = sum(bren_x(:)) + sum(bren_y(:)) + sum(bren_z(:));
                
        end

    case 'Tenenbaum'
        
        % compute sharpness metric based on the Tenenbaum gradient
        switch numDim(im)
            case 2
                
                % define the 2D sobel gradient operator
                sobel = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
                
                % compute metric
                s = conv2(sobel, im).^2 + conv2(sobel.', im).^2;
                s = sum(s(:)); 
                
            case 3
                
                % define the 3D sobel gradient operator
                sobel3D(:, :, 1) = [1, 2, 1; 2, 4, 2; 1, 2, 1];
                sobel3D(:, :, 2) = zeros(3);
                sobel3D(:, :, 3) = -sobel3D(:, :, 1);
                
                % compute metric
                s = convn(im, sobel3D).^2 + ...
                    convn(im, permute(sobel3D, [3, 1, 2])).^2 + ...
                    convn(im, permute(sobel3D, [2, 3, 1])).^2;
                s = sum(s(:));
                
        end        
        
    case 'NormVariance'
        
        % compute sharpness metric based on the normalised variance
        mu = mean(im(:));
        s = sum((im(:) - mu).^2) ./ mu;
        
end          