function test_pass = createCWSignals_extractAmpPhase(plot_comparisons, ~)
% DESCRIPTION:
%       Unit test to test the functions createCWSignals and
%       extractAmpPhase.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 11th January 2017
%       last update - 28th April 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby

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

% signal to noise used for noisy test
SNR = 30;

% set comparison threshold
EXACT_COMPARISON_THRESH = 1e-13;
NOISY_COMPARISON_THRESH = 0.2;

% =========================================================================
% TEST 1: EXACT, 2D
% =========================================================================

% define sampling parameters, using an exact number of samples to make it
% periodic
f = 5e6;
PPP = 10;
dt_periodoc = 1/(PPP*f);
t_array_periodic = (0:(PPP*10-1))*dt_periodoc;

% define amplitude and phase
amp = getWin(9, 'Gaussian');
phase = linspace(0, 2*pi, 9).';

% create 2D signals from the amplitude and phase, without an upramp
cw_signal_2D = createCWSignals(t_array_periodic, f, amp, phase, 0);

% re-extract the amplitude and phase without using a Window, and no zero
% padding
[extracted_amp, extracted_phase] = extractAmpPhase(cw_signal_2D, 1/dt_periodoc, f, 'Window', 'Rectangular', 'FFTPadding', 1);

% unwrap the phase
extracted_phase = unwrap(extracted_phase);

% normalise phase to the first value
extracted_phase = extracted_phase - extracted_phase(1);

% check sizes are the same
if any(size(extracted_amp) ~= size(amp)) || any(size(extracted_phase) ~= size(amp))
    test_pass = false;
end

% check values are the same
if (max(abs(extracted_amp(:) - amp(:))) > EXACT_COMPARISON_THRESH) || ...
   (max(abs(extracted_phase(:) - phase(:))) > EXACT_COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    % plot signals
    figure;
    stackedPlot(cw_signal_2D);
    
    % plot amplitude
    figure;
    subplot(2, 2, 1);
    plot(amp, 'b-');
    hold on;
    plot(extracted_amp, 'r--.');
    title('Amplitudes');
    legend('Original', 'Extracted', 'Location', 'Best');
    subplot(2, 2, 3);
    plot(amp - extracted_amp, 'k-');
    title('Difference');
    
    % plot phase
    subplot(2, 2, 2);
    plot(phase, 'b-');
    hold on;
    plot(extracted_phase, 'r--.');
    title('Amplitudes');
    legend('Original', 'Extracted', 'Location', 'Best');
    subplot(2, 2, 4);
    plot(phase - extracted_phase, 'k-');
    title('Difference');
    
end

% =========================================================================
% TEST 2: 1D DIFFERENT ORIENTATIONS
% =========================================================================

% define sampling parameters
f = 5e6;
Fs = 63e6;
t_array = (0:1/Fs:10*1/f);

% define amplitude and phase
amp = 1;
phase = 0;

% create 1D signals from the amplitude and phase
cw_signal_1D = createCWSignals(t_array, f, amp, phase);

% extract the amplitude and phase with default parameters with the data in
% different directions
[extracted_amp_dim1, extracted_phase_dim1] = extractAmpPhase(reshape(cw_signal_1D, [], 1, 1, 1), Fs, f, 'Dim', 1);
[extracted_amp_dim2, extracted_phase_dim2] = extractAmpPhase(reshape(cw_signal_1D, 1, [], 1, 1), Fs, f, 'Dim', 2);
[extracted_amp_dim3, extracted_phase_dim3] = extractAmpPhase(reshape(cw_signal_1D, 1, 1, [], 1), Fs, f, 'Dim', 3);
[extracted_amp_dim4, extracted_phase_dim4] = extractAmpPhase(reshape(cw_signal_1D, 1, 1, 1, []), Fs, f, 'Dim', 4);

% check the values
err_amp_12 = abs(extracted_amp_dim1 - extracted_amp_dim2);
err_amp_13 = abs(extracted_amp_dim1 - extracted_amp_dim3);
err_amp_14 = abs(extracted_amp_dim1 - extracted_amp_dim4);
err_phs_12 = abs(extracted_phase_dim1 - extracted_phase_dim2);
err_phs_13 = abs(extracted_phase_dim1 - extracted_phase_dim3);
err_phs_14 = abs(extracted_phase_dim1 - extracted_phase_dim4);

% check values are the same
if (err_amp_12 > EXACT_COMPARISON_THRESH) || ...
        (err_amp_13 > EXACT_COMPARISON_THRESH) || ...
        (err_amp_14 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_12 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_13 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_14 > EXACT_COMPARISON_THRESH)
    test_pass = false;
end

% =========================================================================
% TEST 3: EXACT, 3D
% =========================================================================

% define amplitude and phase
amp = getWin([9, 9], 'Gaussian');
phase = repmat(linspace(0, 2*pi, 9).', [1, 9]);

% create 3D signals from the amplitude and phase, without an upramp
cw_signal_3D = createCWSignals(t_array_periodic, f, amp, phase, 0);

% re-extract the amplitude and phase without using a Window, and no zero
% padding
[extracted_amp, extracted_phase] = extractAmpPhase(cw_signal_3D, 1/dt_periodoc, f, 'Window', 'Rectangular', 'FFTPadding', 1);

% unwrap the phase
extracted_phase = unwrap(extracted_phase);

% normalise phase to the first value
extracted_phase = extracted_phase - extracted_phase(1);

% check sizes are the same
if any(size(extracted_amp) ~= size(amp)) || any(size(extracted_phase) ~= size(amp))
    test_pass = false;
end

% check values are the same
if (max(abs(extracted_amp(:) - amp(:))) > EXACT_COMPARISON_THRESH) || ...
   (max(abs(extracted_phase(:) - phase(:))) > EXACT_COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    % plot amplitude
    figure;
    subplot(2, 3, 1);
    imagesc(amp);
    colorbar;
    axis image;
    title('Original Amplitude');
    subplot(2, 3, 2);
    imagesc(extracted_amp);
    colorbar;
    axis image;
    title('Extracted Amplitude');
    subplot(2, 3, 3);
    imagesc(amp - extracted_amp);
    colorbar;
    axis image;
    title('Difference');
    
    % plot phase
    subplot(2, 3, 4);
    imagesc(phase);
    colorbar;
    axis image;
    title('Original Phase');
    subplot(2, 3, 5);
    imagesc(extracted_phase);
    colorbar;
    axis image;
    title('Extracted Phase');
    subplot(2, 3, 6);
    imagesc(phase - extracted_phase);
    colorbar;
    axis image;
    title('Difference');
    
end

% =========================================================================
% TEST 4: NOISY, 2D
% =========================================================================

% define sampling parameters
f = 5e6;
Fs = 63e6;
t_array = (0:1/Fs:10*1/f);

% define amplitude and phase
amp = getWin(9, 'Gaussian');
phase = linspace(0, 2*pi, 9).';

% create 2D signals from the amplitude and phase, without an upramp
cw_signal_2D = createCWSignals(t_array, f, amp, phase, 0);

% add some noise
cw_signal_2D = addNoise(cw_signal_2D, SNR);

% re-extract the amplitude and phase with default parameters
[extracted_amp, extracted_phase] = extractAmpPhase(cw_signal_2D, Fs, f);

% unwrap the phase
extracted_phase = unwrap(extracted_phase);

% normalise phase to the middle value (largest amplitude)
phase = phase - phase(5);
extracted_phase = extracted_phase - extracted_phase(5);

% check sizes are the same
if any(size(extracted_amp) ~= size(amp)) || any(size(extracted_phase) ~= size(amp))
    test_pass = false;
end

% check values are the same
if (max(abs(extracted_amp(:) - amp(:))) > NOISY_COMPARISON_THRESH) || ...
   (max(abs(extracted_phase(:) - phase(:))) > NOISY_COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    % plot signals
    figure;
    stackedPlot(cw_signal_2D);
    
    % plot amplitude
    figure;
    subplot(2, 2, 1);
    plot(amp, 'b-');
    hold on;
    plot(extracted_amp, 'r--.');
    title('Amplitudes');
    legend('Original', 'Extracted', 'Location', 'Best');
    subplot(2, 2, 3);
    plot(amp - extracted_amp, 'k-');
    title('Difference');
    
    % plot phase
    subplot(2, 2, 2);
    plot(phase, 'b-');
    hold on;
    plot(extracted_phase, 'r--.');
    title('Amplitudes');
    legend('Original', 'Extracted', 'Location', 'Best');
    subplot(2, 2, 4);
    plot(phase - extracted_phase, 'k-');
    title('Difference');
    
end

% =========================================================================
% TEST 5: NOISY, 3D
% =========================================================================

% define amplitude and phase
amp = getWin([9, 9], 'Gaussian');
phase = repmat(linspace(0, 2*pi, 9).', [1, 9]);

% create 3D signals from the amplitude and phase, without an upramp
cw_signal_3D = createCWSignals(t_array, f, amp, phase, 0);

% add some noise
cw_signal_3D = addNoise(cw_signal_3D, SNR);

% re-extract the amplitude and phase without using a Window, and no zero
% padding
[extracted_amp, extracted_phase] = extractAmpPhase(cw_signal_3D, Fs, f);

% unwrap the phase
extracted_phase = unwrap(extracted_phase);

% normalise phase to the middle value (largest amplitude)
phase = phase - phase(5, 5);
extracted_phase = extracted_phase - extracted_phase(5, 5);

% check sizes are the same
if any(size(extracted_amp) ~= size(amp)) || any(size(extracted_phase) ~= size(amp))
    test_pass = false;
end

% check values are the same
if (max(abs(extracted_amp(:) - amp(:))) > NOISY_COMPARISON_THRESH) || ...
   (max(abs(extracted_phase(:) - phase(:))) > NOISY_COMPARISON_THRESH)
    test_pass = false;
end

% plot
if plot_comparisons
    
    % plot amplitude
    figure;
    subplot(2, 3, 1);
    imagesc(amp);
    colorbar;
    axis image;
    title('Original Amplitude');
    subplot(2, 3, 2);
    imagesc(extracted_amp);
    colorbar;
    axis image;
    title('Extracted Amplitude');
    subplot(2, 3, 3);
    imagesc(amp - extracted_amp);
    colorbar;
    axis image;
    title('Difference');
    
    % plot phase
    subplot(2, 3, 4);
    imagesc(phase);
    colorbar;
    axis image;
    title('Original Phase');
    subplot(2, 3, 5);
    imagesc(extracted_phase);
    colorbar;
    axis image;
    title('Extracted Phase');
    subplot(2, 3, 6);
    imagesc(phase - extracted_phase);
    colorbar;
    axis image;
    title('Difference');
    
end

% =========================================================================
% TEST 6: 4D DIFFERENT ORIENTATIONS
% =========================================================================

% define sampling parameters
f = 5e6;
Fs = 63e6;
t_array = (0:1/Fs:10*1/f);

% define amplitude and phase
amp = rand(5, 5, 5);
phase = rand(5, 5, 5);

% create 1D signals from the amplitude and phase
cw_signal_4D = createCWSignals(t_array, f, amp, phase);

% extract the amplitude and phase with default parameters with the data in
% different directions
[extracted_amp_dim1, extracted_phase_dim1] = extractAmpPhase(permute(cw_signal_4D, [1 2 3 4]), Fs, f, 'Dim', 1);
[extracted_amp_dim2, extracted_phase_dim2] = extractAmpPhase(permute(cw_signal_4D, [4 1 2 3]), Fs, f, 'Dim', 2);
[extracted_amp_dim3, extracted_phase_dim3] = extractAmpPhase(permute(cw_signal_4D, [3 4 1 2]), Fs, f, 'Dim', 3);
[extracted_amp_dim4, extracted_phase_dim4] = extractAmpPhase(permute(cw_signal_4D, [2 3 4 1]), Fs, f, 'Dim', 4);

% put data back in correct direction
extracted_amp_dim2 = ipermute(extracted_amp_dim2, [4 1 2 3]);
extracted_amp_dim3 = ipermute(extracted_amp_dim3, [3 4 1 2]);
extracted_amp_dim4 = ipermute(extracted_amp_dim4, [2 3 4 1]);
extracted_phase_dim2 = ipermute(extracted_phase_dim2, [4 1 2 3]);
extracted_phase_dim3 = ipermute(extracted_phase_dim3, [3 4 1 2]);
extracted_phase_dim4 = ipermute(extracted_phase_dim4, [2 3 4 1]);

% check the values
err_amp_12 = max(abs(extracted_amp_dim1(:) - extracted_amp_dim2(:)));
err_amp_13 = max(abs(extracted_amp_dim1(:) - extracted_amp_dim3(:)));
err_amp_14 = max(abs(extracted_amp_dim1(:) - extracted_amp_dim4(:)));
err_phs_12 = max(abs(extracted_phase_dim1(:) - extracted_phase_dim2(:)));
err_phs_13 = max(abs(extracted_phase_dim1(:) - extracted_phase_dim3(:)));
err_phs_14 = max(abs(extracted_phase_dim1(:) - extracted_phase_dim4(:)));

% check values are the same
if (err_amp_12 > EXACT_COMPARISON_THRESH) || ...
        (err_amp_13 > EXACT_COMPARISON_THRESH) || ...
        (err_amp_14 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_12 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_13 > EXACT_COMPARISON_THRESH) || ...
        (err_phs_14 > EXACT_COMPARISON_THRESH)
    test_pass = false;
end