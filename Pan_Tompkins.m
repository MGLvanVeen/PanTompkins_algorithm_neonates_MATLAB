function [peaks_idx, peaks_data, THRESHOLD_I1, smoothed_data] = Pan_Tompkins(fs, data, ecg_column, k, plot_flag)
% Pan_Tompkins
% 入力:
%   fs         - サンプリング周波数 (Hz)
%   data       - 生のデータ行列
%   ecg_column - データ内のECG信号が含まれる列番号
%   k          - 閾値の何倍を超えたピークを閾値更新から除外する条件
%   plot_flag  - 図を生成する場合は1, 生成しない場合は0
% 出力:
%   peaks_idx - 検出されたR波のインデックス
%   peaks_data - 検出されたR波の振幅(大きさ)
%   smoothed_data - 前処理されたECG信号
%   THRESHOLD_I1 - 移動閾値信号(第一閾値，第二閾値はこれの半分)
% ブロック化を行いたいときにはpeaks_idxが出力されるので求めたい時間(fs*time)以上以下で求められる

%% データの前処理
data = data(:, ecg_column); % 指定列からECGデータを抽出
data = data(isfinite(data)); % NaN値を削除
data = data(:);
% バンドパスフィルタの適用
[b, a] = butter(5, [5 / (fs / 2), 15 / (fs / 2)], 'bandpass'); % フィルタ係数
filtered_data = filtfilt(b, a, data);
% 微分フィルタ
derivative_filter = [1 2 0 -2 -1] .* (1 / 8) * fs;
differentiated_data = filtfilt(derivative_filter, 1, filtered_data);
% 平方と移動平均
squared_data = differentiated_data .^ 2;
smoothed_data = movmean(squared_data, round(0.15 * fs));

%% 初期閾値設定
refractory_period=round(0.2*fs);%不応期
initial_segment = smoothed_data(1:2 * fs);%学習フェーズ
[pks, ~] = findpeaks(initial_segment, 'MinPeakDistance', refractory_period);
initial_SPKI = mean(pks);
initial_NPKI = mean(pks(pks < prctile(pks, 50)));

%% ピーク検出と閾値設定
% ピーク検出とリアルタイム再検索
signal_peaks_idx = [];
SPKI = initial_SPKI;
NPKI = initial_NPKI;
retrieved_peaks_idx = []; % 再検索されたR波ピークを格納する配列
signal_pks_data = [];
t_wave_detected = [];
t_wave_amplitudes = [];
THRESHOLD_I1 = zeros(1, length(smoothed_data)); % 閾値用の信号
threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
threshold_I2 = 0.5 * threshold_I1;

% ピーク検出(R波候補)
[pks, locs] = findpeaks(smoothed_data, 'MinPeakDistance', refractory_period); 

% インデックスごとのリアルタイム再検索とT波識別を行う閾値更新
for i = 1:length(locs)
    % 再検索の実行
    if ~isempty(signal_peaks_idx)
        % RR間隔を計算
        RR_interval = diff(signal_peaks_idx);
        % RRaverage1: 規則的な心拍リズムを考慮した平均RR間隔
        RRaverage1 = mean(RR_interval(max(1, end-7):end));
        % RRaverage2: 限定範囲内 (92%～116%) の平均RR間隔
        valid_RR = RR_interval(RR_interval >= 0.92 * RRaverage1 & RR_interval <= 1.16 * RRaverage1);
        if ~isempty(valid_RR)
            RRaverage2 = mean(valid_RR);
        else
            RRaverage2 = RRaverage1; % 規則的な平均を代用
        end
        % 見逃しチェック: 166%以内にR波がない場合、再検索
        if (locs(i) - signal_peaks_idx(end)) > 1.66 * RRaverage2
            search_start = signal_peaks_idx(end) + refractory_period;
            search_end = locs(i) - refractory_period;
            range_length = search_end - search_start + 1; % 再検索範囲の長さを計算
            re_locs=[];
            if range_length > refractory_period
                % 再検索範囲が不応期より大きい場合
                [re_pks, re_locs] = findpeaks(smoothed_data(search_start:search_end),  'MinPeakDistance', refractory_period);
            end
            % 再検索範囲でピークを検出
            if ~isempty(re_locs)
                % 第二閾値以上のピークを判定
                above_threshold_idx = find(re_pks >= threshold_I2);
                for j = above_threshold_idx'
                    current_peak = re_pks(j);
                    current_loc = re_locs(j) + search_start ;
                    % 再検索されたピークを追加
                    retrieved_peaks_idx = [retrieved_peaks_idx; current_loc];
                    % ピークを追加し、閾値を更新
                    signal_peaks_idx = [signal_peaks_idx; current_loc];
                    signal_pks_data = [signal_pks_data; current_peak];
                    SPKI = 0.125 * current_peak + 0.875 * SPKI;
                    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
                    threshold_I2 = 0.5 * threshold_I1;
                end
            end
        end
    end
    
    % T波の判定フラグ
    is_t_wave = false;
   % T波識別 (現在のピークと直前のピークの比較)
    if ~isempty(signal_peaks_idx)
        last_r_peak = signal_peaks_idx(end); % 最後に検出されたR波ピーク
        if (locs(i) - last_r_peak) < round(0.36 * fs)
            % 勾配を比較してT波判定
            window_size = round(0.075 * fs); % 75msの窓幅
            % 前のピークの勾配を計算 (範囲外を防止)
            if last_r_peak > window_size
                prev_slope = max(diff(smoothed_data(last_r_peak-window_size:last_r_peak)));
            else
                prev_slope = max(diff(smoothed_data(1:last_r_peak))); % データ先頭から現在位置まで使用
            end
            % 現在のピークの勾配を計算 (範囲外を防止)
            if locs(i) > window_size
                curr_slope = max(diff(smoothed_data(locs(i)-window_size:locs(i))));
            else
                curr_slope = max(diff(smoothed_data(1:locs(i)))); % データ先頭から現在位置まで使用
            end
            % T波判定
            if curr_slope < 0.5 * prev_slope
                t_wave_detected = [t_wave_detected; locs(i)]; % T波として識別されたピーク
                t_wave_amplitudes = [t_wave_amplitudes; pks(i)]; % T波として識別されたピーク
                is_t_wave = true;
            end
        end
    end

    % ピーク分類 (T波の場合はノイズとして扱う)
    if is_t_wave
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    elseif pks(i) > k * threshold_I1
        % k倍を超えたピーク：ピークとして検出するが閾値更新には使用しない
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    elseif pks(i) > threshold_I1
        % 正常なR波ピーク
        SPKI = 0.125 * pks(i) + 0.875 * SPKI;
        signal_peaks_idx = [signal_peaks_idx; locs(i)];
        signal_pks_data = [signal_pks_data; pks(i)];
    else
        % ノイズピーク
        NPKI = 0.125 * pks(i) + 0.875 * NPKI;
    end 
    % 閾値更新
    threshold_I1 = NPKI + 0.25 * (SPKI - NPKI);
    threshold_I2 = 0.5 * threshold_I1;
    % 信号に閾値を反映
    if i == 1
        THRESHOLD_I1(1:locs(i)) = threshold_I1;
    elseif i == length(locs)
        THRESHOLD_I1(locs(i-1)+1:end) = threshold_I1;
    else
        THRESHOLD_I1(locs(i-1)+1:locs(i)) = threshold_I1;
    end
end
% ===== 再検索とT波識別後のデータを保存 =====
peaks_idx = signal_peaks_idx;
peaks_data = smoothed_data(peaks_idx);

%% プロット
if plot_flag == 1
    time = (0:length(smoothed_data) - 1) / fs; % 時間軸作成
    figure;
    plot(time, smoothed_data, 'b', 'DisplayName', 'Smoothed Signal'); % 平滑化された信号
    hold on;

    % 再検出されたR波をプロット
    if ~isempty(retrieved_peaks_idx)
        valid_indices = retrieved_peaks_idx(retrieved_peaks_idx > 0 & retrieved_peaks_idx <= length(smoothed_data)); 
        retrieved_amplitudes = smoothed_data(valid_indices-1);% 時間軸にあわせるためvalid_indices事態が時間軸に対応しているため(0秒からなので)
        plot((valid_indices-2) / fs, retrieved_amplitudes, 'go', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Re-detected)');% -2は時間軸に対応させるための調整　
    end

    % 第一閾値で検出されたR波をプロット
    if ~isempty(peaks_idx)
        initial_peaks = setdiff(peaks_idx, retrieved_peaks_idx); % 再検索ピークを除外
        initial_peaks = initial_peaks(initial_peaks > 0 & initial_peaks <= length(smoothed_data)); 
        if ~isempty(initial_peaks)
            plot((initial_peaks-1) / fs, smoothed_data(initial_peaks), 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'R-peaks (Initial)');%0秒がデータだと1つ目だから-1
        end
    end
    % 第二閾値以上のT波をプロット
    if ~isempty(t_wave_detected)
        % 保存用の配列を初期化
        valid_t_wave_times = [];
        valid_t_wave_amplitudes = [];
        for i = 1:length(t_wave_detected)
            if t_wave_amplitudes(i) > THRESHOLD_I1(t_wave_detected(i)) / 2
                % 有効なT波のインデックスと振幅を保存
                valid_t_wave_times = [valid_t_wave_times, t_wave_detected(i)];
                valid_t_wave_amplitudes = [valid_t_wave_amplitudes, t_wave_amplitudes(i)];
            end
        end
        if ~isempty(valid_t_wave_times)
            plot((valid_t_wave_times-1) / fs, valid_t_wave_amplitudes, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'T-wave (Excluded)'); %0秒がデータだと1つ目だから-1
        end
    end
    % 閾値をプロット
    plot(time, THRESHOLD_I1, 'm--', 'LineWidth', 1, 'DisplayName', 'Threshold I1');
    plot(time, THRESHOLD_I1 / 2, 'g--', 'LineWidth', 1, 'DisplayName', 'Threshold I2');
    title('ECG Signal with Detected R-peaks and T-wave Exclusion');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('show', 'Location', 'best');
    grid on;
    hold off;
end
end