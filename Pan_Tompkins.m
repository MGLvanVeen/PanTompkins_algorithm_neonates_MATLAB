function [peaks_idx, peaks_data, THRESHOLD_I1, smoothed_data, clean_time] = Pan_Tompkins(fs, data, time, k, plot_flag) 
% Pan_Tompkins

% Adapted by MGL van Veen from: "Pan-Tompkins Algorithm Implementation in MATLAB" (https://github.com/Tatsuo5/PanTompkins_algorithm_MATLAB) by Tatsuo Hirata (MIT License)
% This fork (https://github.com/MGLvanVeen/PanTompkins_algorithm_neonates_MATLAB) contains minor adaptations for dataset/time handling and robustness in re-detection and T-wave logic.
% Last update: 27-10-2025

% SPDX-License-Identifier: MIT
% Copyright (c) 2024 Tatsuo Hirata
% Copyright (c) 2025 M. G. L. van Veen

% Summary of adaptations (see README and CHANGELOG for details):
% - Added explicit time input and clean_time output for NaN-safe alignment.  % ADAPTED
% - Robust re-detection bounds and dynamic MinPeakDistance.                  % ADAPTED
% - Adaptive T-wave window based on recent RR average (with fallback).       % ADAPTED
% - Interface simplified: pass ECG vector directly (no ecg_column).          % ADAPTED

% Input:
%   fs         - Sampling frequency (Hz)
%   data       - Raw data vector or matrix (ECG channel already selected) % ADAPTED
%   time       - Time vector aligned with data (same length) % ADAPTED
%   k          - Condition to exclude peaks from threshold update if they exceed k times the threshold
%   plot_flag  - 1 to generate plots, 0 to not generate plots
% Output:
%   peaks_idx       - Detected R-wave indices
%   peaks_data      - Amplitudes (magnitudes) of the detected R-waves
%   smoothed_data   - Preprocessed ECG signal
%   THRESHOLD_I1    - Moving threshold signal (first threshold, with the second threshold being half of this)
%   clean_time      - Time vector after NaN removal (aligned with cleaned data) % ADAPTED
% When block processing is desired, peaks_idx is output so that R-waves can be determined for a given time (fs*time) interval

%% Data Preprocessing
valid_idx = isfinite(data);           % logical index of valid entries %ADAPTED build robust mask
data = data(valid_idx);               % clean signal by removing Nan values %ADAPTED
clean_time = time(valid_idx);         % clean time vector %ADAPTED use mask so the time vector stays aligned with the cleaned signal
data = data(:);                       % ADAPTED ensure column vector for consistent filtering downstream
% Apply bandpass filter
[b, a] = butter(5, [5 / (fs / 2), 15 / (fs / 2)], 'bandpass'); % Filter coefficients
filtered_data = filtfilt(b, a, data);
% Differentiation filter
derivative_filter = [1 2 0 -2 -1] .* (1 / 8) * fs;
differentiated_data = filtfilt(derivative_filter, 1, filtered_data);
% Squaring and moving average
squared_data = differentiated_data .^ 2;
smoothed_data = movmean(squared_data, round(0.15 * fs));

%% Initial Threshold Setting
refractory_period = round(0.2 * fs); % Refractory period
initial_segment = smoothed_data(1:2 * fs); % Learning phase
[pks, ~] = findpeaks(initial_segment, 'MinPeakDistance', refractory_period);
initial_SPKI = mean(pks);
initial_NPKI = mean(pks(pks < prctile(pks, 50)));

%% Peak Detection and Threshold Setting
% Peak detection and real-time re-detection
signal_peaks_idx = [];
SPKI = initial_SPKI;
NPKI = initial_NPKI;
retrieved_peaks_idx = []; % Array to store re-detected R-wave peaks
signal_pks_data = [];
t_wave_detected = [];
t_wave_amplitudes = [];
THRESHOLD_I1 = zeros(1, length(smoothed_data)); % Signal for threshold
threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
threshold_I2 = 0.5 * threshold_I1;
% Peak detection (R-wave candidates)
[pks, locs] = findpeaks(smoothed_data, 'MinPeakDistance', refractory_period); 
% Real-time re-detection and T-wave identification with threshold update for each index
for i = 1:length(locs)
    % Execute re-detection
    if ~isempty(signal_peaks_idx)
        % Calculate RR intervals
        RR_interval = diff(signal_peaks_idx);
        % RRaverage1: Average RR interval considering regular heart rhythm
        RRaverage1 = mean(RR_interval(max(1, end-7):end));
        % RRaverage2: Average RR interval within the range (92% to 116%)
        valid_RR = RR_interval(RR_interval >= 0.92 * RRaverage1 & RR_interval <= 1.16 * RRaverage1);
        if ~isempty(valid_RR)
            RRaverage2 = mean(valid_RR);
        else
            RRaverage2 = RRaverage1; % Substitute with the regular average
        end
        % Missed beat check: if no R-wave within 166%, re-detect
        if (locs(i) - signal_peaks_idx(end)) > 1.66 * RRaverage2
            search_start = signal_peaks_idx(end) + refractory_period;
            search_end = locs(i) - refractory_period;
            re_locs = []; % ADAPTED
            if search_end > search_start % ADAPTED guard against empty/negative ranges
                range = smoothed_data(search_start:search_end); % ADAPTED slice once to avoid repeated indexing
                x_span = search_end - search_start; % ADAPTED span for constraints
                if length(range) > 2 && x_span > 1 % ADAPTED ensure enough samples to detect peaks and the re-detection range is larger than the refractory period
                    mpd = min(refractory_period, x_span - 1); % ADAPTED cap MinPeakDistance to range length to prevent errors
                    [re_pks, re_locs] = findpeaks(range, 'MinPeakDistance', mpd); % ADAPTED robust re-search even in short gaps
                end
            end
            % Detect peaks in the re-detection range
            if ~isempty(re_locs)
                % Determine peaks above the second threshold
                above_threshold_idx = find(re_pks >= threshold_I2);
                % If there are peaks above the second threshold, process them all
                for j = above_threshold_idx'
                    current_peak = re_pks(j);
                    current_loc = re_locs(j) + search_start;
                    is_t_wave = false;
                    if ~isempty(signal_peaks_idx)
                        last_r_peak = signal_peaks_idx(end); 
                        if (current_loc - last_r_peak) < round(0.36 * fs)
                            window_size = round(0.075 * fs);
                            if last_r_peak > window_size
                                prev_slope = max(diff(smoothed_data(last_r_peak-window_size:last_r_peak)));
                            else
                                prev_slope = max(diff(smoothed_data(1:last_r_peak)));
                            end
                            if current_loc > window_size
                                curr_slope = max(diff(smoothed_data(current_loc-window_size:current_loc)));
                            else
                                curr_slope = max(diff(smoothed_data(1:current_loc)));
                            end
                            if curr_slope < 0.5 * prev_slope
                                t_wave_detected = [t_wave_detected; current_loc];
                                t_wave_amplitudes = [t_wave_amplitudes; current_peak];
                                is_t_wave = true;
                            end
                        end
                    end
                    if is_t_wave
                        NPKI = 0.125 * current_peak + 0.875 * NPKI;
                    else
                        retrieved_peaks_idx = [retrieved_peaks_idx; current_loc];
                        signal_peaks_idx = [signal_peaks_idx; current_loc];
                        signal_pks_data = [signal_pks_data; current_peak];
                        SPKI = 0.125 * current_peak + 0.875 * SPKI;
                    end
                    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
                    threshold_I2 = 0.5 * threshold_I1;
                end
            end
        end
    end
    % T-wave flag
    is_t_wave = false;
    % T-wave identification (comparing the current peak with the previous peak)
    if ~isempty(signal_peaks_idx)
        last_r_peak = signal_peaks_idx(end); % Last detected R-wave peak
        RR_interval = diff(signal_peaks_idx); % ADAPTED compute RR intervals for adaptive timing window
        if ~isempty(RR_interval) % ADAPTED if we have history, use recent mean RR
            RRaverage = mean(RR_interval(max(1, end-7):end)); % ADAPTED average over last up-to-8 beats for stability
        else
            RRaverage = fs * 0.8;  % ADAPTED fallback = 800 ms at given fs to avoid edge cases
        end
        
        if (locs(i) - last_r_peak) < 0.35 * RRaverage % ADAPTED to adaptive RR (heart rate aware) from fixed round (0.36*fs)
            % Compare slopes to determine T-wave
            window_size = round(0.075 * fs); % 75ms window
            % Calculate the slope of the previous peak (preventing out-of-bound errors)
            if last_r_peak > window_size
                prev_slope = max(diff(smoothed_data(last_r_peak-window_size:last_r_peak)));
            else
                prev_slope = max(diff(smoothed_data(1:last_r_peak))); % Use from the start of data to current position
            end
            % Calculate the slope of the current peak (preventing out-of-bound errors)
            if locs(i) > window_size
                curr_slope = max(diff(smoothed_data(locs(i)-window_size:locs(i))));
            else
                curr_slope = max(diff(smoothed_data(1:locs(i)))); % Use from the start of data to current position
            end
            % T-wave determination
            if curr_slope < 0.5 * prev_slope
                t_wave_detected = [t_wave_detected; locs(i)]; % Peak identified as T-wave
                t_wave_amplitudes = [t_wave_amplitudes; pks(i)]; % Amplitude of the peak identified as T-wave
                is_t_wave = true;
            end
        end
    end
    % Peak classification (if T-wave, treat as noise)
    if is_t_wave
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    elseif pks(i) > k * threshold_I1
        % Peak exceeding k times: detected as peak but not used for threshold update
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    elseif pks(i) > threshold_I1
        % Normal R-wave peak
        SPKI = 0.125 * pks(i) + 0.875 * SPKI;
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    else
        % Noise peak
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    end 
    % Update thresholds
    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
    threshold_I2 = 0.5 * threshold_I1;
    % Reflect threshold on signal
    if i == 1
        THRESHOLD_I1(1:locs(i)) = threshold_I1;
    elseif i == length(locs)
        THRESHOLD_I1(locs(i-1)+1:end) = threshold_I1;
    else
        THRESHOLD_I1(locs(i-1)+1:locs(i)) = threshold_I1;
    end
end
% Save data after re-detection and T-wave identification
peaks_idx = signal_peaks_idx;
peaks_data = smoothed_data(peaks_idx);

%% Plot
if plot_flag == 1
    time = (0:length(smoothed_data) - 1) / fs; % ADAPTED regenerate a uniform time axis for plotting
    figure;
    plot(time, smoothed_data, 'b', 'DisplayName', 'Smoothed Signal'); % Plot smoothed signal % ADAPTED plot against uniform time axis
    hold on;
    % Plot re-detected R-waves
    if ~isempty(retrieved_peaks_idx)
        valid_indices = retrieved_peaks_idx(retrieved_peaks_idx > 0 & retrieved_peaks_idx <= length(smoothed_data)); 
        retrieved_amplitudes = smoothed_data(valid_indices-1); % Adjusted for time axis (since indexing starts at 0)
        plot((valid_indices-2) / fs, retrieved_amplitudes, 'go', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Re-detected)'); % ADAPTED plot against uniform time axis
    end
    % Plot R-waves detected by the first threshold
    if ~isempty(peaks_idx)
        initial_peaks = setdiff(peaks_idx, retrieved_peaks_idx); % Exclude re-detected peaks
        initial_peaks = initial_peaks(initial_peaks > 0 & initial_peaks <= length(smoothed_data)); 
        if ~isempty(initial_peaks)
            plot((initial_peaks-1) / fs, smoothed_data(initial_peaks), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Initial)'); % ADAPTED plot against uniform time axis
        end
    end
    % Plot T-waves above the second threshold
    if ~isempty(t_wave_detected)
        valid_t_wave_times = [];
        valid_t_wave_amplitudes = [];
        for i = 1:length(t_wave_detected) 
            if t_wave_amplitudes(i) > THRESHOLD_I1(t_wave_detected(i)) / 2
                valid_t_wave_times = [valid_t_wave_times, t_wave_detected(i)];
                valid_t_wave_amplitudes = [valid_t_wave_amplitudes, t_wave_amplitudes(i)];
            end
        end
        if ~isempty(valid_t_wave_times)
            plot((valid_t_wave_times-1) / fs, valid_t_wave_amplitudes, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'T-wave (Excluded)'); % ADAPTED plot against uniform time axis
        end
    end
    % Plot thresholds
    plot(time, THRESHOLD_I1, 'm--', 'LineWidth', 1, 'DisplayName', 'Threshold I1'); % ADAPTED plot against uniform time axis
    plot(time, THRESHOLD_I1 / 2, 'g--', 'LineWidth', 1, 'DisplayName', 'Threshold I2'); % ADAPTED plot against uniform time axis
    title('ECG Signal with Detected R-peaks and T-wave Exclusion');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
end
end
