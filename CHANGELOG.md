# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2025-10-27
### Added
- New input `time` and new output `clean_time` to maintain alignment after NaN removal.
- NaN/Inf filtering with shared mask applied to both `data` and `time`.
- Adaptive T-wave window using `0.35 × RRaverage` with fallback when insufficient history exists.
- Robust re-detection range handling (bounds checks and dynamic `MinPeakDistance`).
- Inline comments tagged `% ADAPTED` for every change from upstream.
- README notes that reproduce original wording and mark adaptations clearly.
- MIT license header guidance for `.m` files.

### Changed
- Interface no longer requires `ecg_column`; `data` is passed as the ECG vector directly.
- Plotting defaults to a uniform time axis; notes provided to plot against `clean_time` if preferred.
- Minor clarifications and safety checks in re-search and T-wave logic to better handle neonatal ECG.

### Removed
- Original sample dataset and demo code not used in this adaptation.

[1.0.0]: https://github.com/MGLvanVeen/PanTompkins_algorithm_neonates_MATLAB/releases/tag/v1.0.0
