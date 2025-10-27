# Pan-Tompkins Algorithm Implementation in MATLAB (Adapted Version for Neonatal ECG)

**Adapted by:** M. G. L. van Veen  
**Date:** 27 October 2025  
**Original repository:** [Tatsuo Hirata – PanTompkins_algorithm_MATLAB](https://github.com/Tatsuo5/PanTompkins_algorithm_MATLAB)  
**This fork:** [MGL van Veen – PanTompkins_algorithm_neonates_MATLAB](https://github.com/MGLvanVeen/PanTompkins_algorithm_neonates_MATLAB)  
**License:** MIT License (inherits from the original project)

  *This README reproduces the original text and structure of the **Features**, **Usage** and **Algorithm Details** with addition of the changes that were made.*

## Overview 

This repository is a **research fork** of Tatsuo Hirata’s MATLAB implementation of the **Pan-Tompkins algorithm for QRS detection in ECG signals**.  
The algorithm was adapted for use in the analysis of **neonatal ECG data**, described in the publication:
*Cerebral critical closing pressure during neonatal cardiac surgery: a comparison of the 2-point, first harmonic and impedance method.* (2025), by M. G. L. van Veen *et all.* (currently under submission).

The core algorithm logic remains unchanged.  
The modifications focus on:

- Handling of data vectors that include missing (NaN) samples.  
- Incorporation of a corresponding **time vector** for precise alignment.  
- Minor robustness and thresholding adjustments to better handle neonatal ECG morphology.

In the code, adapted lines are indicated by the comment `% ADAPTED`.

## Summary of Adaptations

| Area | Change | Purpose |
|------|---------|----------|
| **Input and Output** | Added input `time` and output `clean_time`. | Align cleaned data with its time axis. |
| **Pre-processing** | Replaced column selection (`ecg_column`) with direct signal input (`data`) and added a finite-sample mask. | Robust cleaning of NaN/Inf values in neonatal recordings. |
| **Re-detection step** | Added checks for short intervals and adaptive `MinPeakDistance`. | Prevent errors on shorter RR intervals. |
| **T-wave identification** | Made window adaptive to current RR-average (`0.35 × RRavg` instead of fixed `0.36 × fs`). | Handle faster neonatal heart rates and reduce false T-wave labelling. |
| **Plotting** | Option to plot against uniform or cleaned time axis. | Maintain alignment between processed signal and timestamps. |
| **Sample data** | Removed `sample_ecg.csv` and demo code. | Not used or validated for this adaptation. |


## Features

- Detects QRS complexes in ECG signals using the Pan-Tompkins algorithm.
- Includes preprocessing steps such as bandpass filtering, differentiation, squaring, and moving average smoothing.
- Executes re-search and T-wave identification as needed, and performs adaptive threshold updates after each peak detection.
- Excludes abnormal peaks from threshold updates based on user-defined criteria.
- Accurately identifies QRS complexes by distinguishing R-waves from T-waves.
- Provides customizable visualization of R-peaks, adaptive thresholds, and T-wave exclusion.
- **_Adapted:_** Improved handling of datasets containing NaN values and non-uniform sampling intervals.  
- **_Adapted:_** Modified T-wave detection and re-search logic for higher heart rates in neonatal ECGs.  
- **_Adapted:_** Added support for time vectors (`time`) to maintain alignment between cleaned ECG and time after NaN removal.


## Usage

### Input

- **Sampling frequency (fs)**: The frequency at which the ECG signal is sampled, specified in Hz.
- **Raw data (data)**: A vector containing the ECG signal.
- **Time (time)**: A time vector aligned with the ECG data (same length as `data`).
- **Threshold multiplier (k)**: A user-defined value to determine how abnormal peaks are excluded from threshold updates.  
  Recommended values are approximately 50–100. If you do not wish to exclude abnormal peaks from threshold updates, set `k` to a very high value (e.g., around 1 million).  
  This parameter specifies that when a peak's amplitude exceeds the current threshold by a factor of `k`, it will still be detected as a peak, but its amplitude will not be used to update the adaptive threshold.
- **Plot flag (plot_flag)**: A flag to enable or disable plot output. Set to 1 for visualization or 0 to disable plotting.

### Output

The `Pan_Tompkins` function returns the following outputs:

- **peaks_idx**: Detected R-wave indices.
- **peaks_data**: Amplitudes (magnitudes) of the detected R-waves.
- **smoothed_data**: Preprocessed ECG signal.
- **THRESHOLD_I1**: Moving threshold signal (first threshold, with the second threshold being half of this).
- **clean_time**: Time vector after NaN removal, aligned with the cleaned ECG signal. 


## Algorithm Details

The Pan-Tompkins algorithm processes ECG signals through the following steps:

1. **Preprocessing**:
   - **NaN filtering**: NaN filtering and synchronization of the corresponding time vector before preprocessing.
   - **Bandpass filtering**: Removes noise and baseline wander by retaining frequencies in the range of 5-15 Hz.
   - **Differentiation filter**: Highlights the slope of QRS complexes using the differentiation filter proposed by Pan and Tompkins:

     $$
     H(z) = \frac{1}{8} \left( -z^{-2} - 2z^{-1} + 2z^{1} + z^{2} \right)
     $$

   - **Squaring**: Amplifies high-frequency components and suppresses negative amplitudes to emphasize QRS peaks.
   - **Moving average smoothing**: Extracts the signal envelope to further highlight QRS peaks and reduce noise.

2. **Peak Detection**:
   - Identifies R-peaks using adaptive thresholds.
   - Handling of short RR intervals and inclusion of safety checks to avoid indexing errors in short neonatal ECG segments.

   **Note:** The initial detected R-peak data may contain false detections due to the adaptation of the threshold. It is recommended to exclude the first detected R-peak(s) when analyzing the results.

3. **Re-search**:
   - If a QRS complex is not detected within 166% of the average RR interval (calculated from the most recent 8 detected R-peaks), the algorithm searches for the largest peak in that interval.
   - The re-search process uses a secondary threshold (THRESHOLD_I2) to increase sensitivity while ensuring false detections are minimized.
   - This step ensures robustness in detecting QRS complexes, even in cases of irregular heartbeats or missed peaks.
   - Adaptive handling of `MinPeakDistance` in the re-search step and improved boundary checks to prevent errors in short signals.

4. **T-Wave Identification**:
   - When an RR interval is less than 360 ms and greater than the 200 ms latency period, the algorithm checks whether the current peak is a T-wave or a QRS complex:
     - The maximum slope of the peak is compared to the maximum slope of the preceding QRS complex.
     - If the maximum slope is less than half of the preceding QRS maximum slope, it is classified as a T-wave.
     - Otherwise, it is classified as a QRS complex.
     - Adaptive threshold (0.35 × average RR interval) to reflect neonatal heart rate variability and prevent T-wave misclassification.

5. **Threshold Updates**:
   - **Initial Threshold Calculation**:
     - During the learning phase (first 2 seconds of the signal), the algorithm calculates the initial thresholds:
       - **SPKI (signal peak indicator)**: The mean of peaks detected in the learning phase.
       - **NPKI (noise peak indicator)**: The mean of peaks below the 50th percentile in the learning phase.
   - **Threshold Update Equations**:
     - After every detected R-peak, the thresholds are updated as follows:

       THRESHOLD\_I1 = NPKI + 0.25 (SPKI - NPKI)

       THRESHOLD\_I2 = 0.5 THRESHOLD\_I1

     - Where:
       - **SPKI** is updated with the amplitude of detected signal peaks.
       - **NPKI** is updated with the amplitude of detected noise peaks.

6. **Enhancements**:
   - Excludes abnormal peaks from threshold updates based on user-defined criteria.
   - Provides customizable visualization of results, including detected R-peaks and adaptive thresholds.
   - Visualization  to use the cleaned `time` vector for accurate plotting after NaN removal.

## Sample Code

In the original repository, sample code and data were included. These were not used for development or testing of this adapted version and thus are not included here.

## License

This project is licensed under the same MIT License as the original project. 

© 2024 Tatsuo Hirata  
© 2025 M. G. L. van Veen  

You are free to use, modify, and distribute this code, provided proper attribution is given. See the LICENSE file in the repository for more details.

## Citation

If you use this code or an adaptation of it, please cite both the original and this fork:

Hirata, T. (2024). Pan-Tompkins Algorithm Implementation in MATLAB.
GitHub repository: https://github.com/Tatsuo5/PanTompkins_algorithm_MATLAB

van Veen, M. G. L. (2025). Adapted Pan-Tompkins Algorithm for Neonatal ECG in MATLAB.
GitHub repository: https://github.com/MGLvanVeen/PanTompkins_algorithm_neonates_MATLAB

## Disclaimer

This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, or non-infringement. The authors are not responsible for any damage or loss caused by the use of this software. Users are encouraged to validate results independently.
