%% Example of Using the Pan-Tompkins Function
%% Step 1: Run the Pan-Tompkins Algorithm
% Load the data
data = readmatrix("sample_data.CSV");
fs = 500; % [Hz]
data=data(3:end,:); % Exclude header information
% Run the Pan-Tompkins algorithm
[peaks_idx, peaks_data, THRESHOLD_I1, smoothed_data] = Pan_Tompkins(fs, data, 2, 100, 1);

%% Step 2: Segmentation into Blocks
R_flag = zeros(1, length(smoothed_data)); % Array for flags (initialized to all zeros)
% Mark the R-peak positions with flags
R_flag(peaks_idx) = 1;

% Time intervals and repetitions
interval1_Time = 15; % Duration of interval 1 [sec]
interval2_Time = 10; % Duration of interval 2 [sec]
interval3_Time = 5; % Duration of interval 3 [sec]
num_repeats = 2; % Number of repetitions
interval1_samples = interval1_Time * fs; % Samples for interval 1
interval2_samples = interval2_Time * fs; % Samples for interval 2
interval3_samples = interval3_Time * fs; % Samples for interval 3

% Total samples per cycle
cycle_samples = interval1_samples + interval2_samples + interval3_samples;

% Create a cell array to store data for interval 2 from each repetition
Block_data = cell(1, num_repeats);

% Extract interval 2 data for each repetition
for repeat = 1:num_repeats
    cycle_start = (repeat - 1) * cycle_samples + 1;
    interval2_start = cycle_start + interval1_samples; 
    interval2_end = interval2_start + interval2_samples - 1; 
    Block_data{1, repeat} = R_flag(interval2_start:interval2_end);
    % Calculate RR intervals within interval 2
    r_peak_indices = find(Block_data{1, repeat} == 1);
    rr_intervals_samples = diff(r_peak_indices);
    focus_rr_intervals{1, repeat} = rr_intervals_samples / fs;
end
