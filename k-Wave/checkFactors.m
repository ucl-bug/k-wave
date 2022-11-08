function checkFactors(min_number, max_number)
%CHECKFACTORS Return the maximum prime factor for a range of numbers.
%
% DESCRIPTION:
%     checkFactors loops through the given range of numbers and finds the
%     numbers with the smallest maximum prime factors. This allows suitable
%     grid sizes to be selected to maximise the speed of the FFT (this is
%     fastest for FFT lengths with small prime factors). The output is
%     printed to the command line, and a plot of the factors is generated.
%    
% USAGE:
%     checkFactors(min_number, max_number)
%
% INPUTS:
%     min_number    - integer specifying the lower bound of values to test
%     max_number    - integer specifying the upper bound of values to test
%
% ABOUT:
%     author        - Bradley Treeby
%     date          - 20th April 2011
%     last update   - 27th February 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% extract factors
facs = zeros(1, max_number - min_number);
fac_max = facs;
for index = min_number:max_number
    facs(index - min_number + 1) = length(factor(index));
    fac_max(index - min_number + 1) = max(factor(index));
end

% plot factors
figure;
subplot(2, 1, 1), bar(min_number:max_number, facs);
set(gca, 'XLim', [(min_number -0.5) (max_number + 0.5)]);
title('number of factors');
subplot(2, 1, 2), bar(min_number:max_number, fac_max);
set(gca, 'XLim', [(min_number -0.5) (max_number + 0.5)]);
title('largest factor');

% compute best factors in range
disp('Numbers with a maximum prime factor of 2');
ind = min_number + find(fac_max == 2) - 1;
disp(num2str(ind));
disp('Numbers with a maximum prime factor of 3');
ind = min_number + find(fac_max == 3) - 1;
disp(num2str(ind));
disp('Numbers with a maximum prime factor of 5');
ind = min_number + find(fac_max == 5) - 1;
disp(num2str(ind));
disp('Numbers with a maximum prime factor of 7');
ind = min_number + find(fac_max == 7) - 1;
disp(num2str(ind));
disp('Numbers to avoid (prime numbers)');
nums = min_number:max_number;
disp(num2str(nums(fac_max == nums)));