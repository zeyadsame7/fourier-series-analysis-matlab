# Fourier Series Analysis in MATLAB ("Signal Symphony")

A MATLAB study of Fourier series: computing coefficients numerically, reconstructing signals, applying Fourier-series properties, and showing aliasing under sampling.

Team project for *ECE 231 – Signals and Systems* at Egypt University of Informatics, by Zeyad Sameh and Ahmed Moataz.

## What it does

1. **Fourier series coefficients and reconstruction.** Computes the complex coefficients `C_k` of a periodic square wave (T = 8 s, 50% duty cycle) by numerical integration, then rebuilds the signal from N = 25 and N = 75 harmonics. The Gibbs phenomenon at the discontinuities is clearly visible.
2. **Differentiation property.** Differentiates the signal in time and checks the Fourier-series property `C_k' = jkω₀·C_k`, plotting the magnitude and phase spectra of the signal and its derivative.
3. **Sampling and aliasing.** Samples `cos(2π·100t) + cos(2π·50t)` at 1000 Hz, 200 Hz and 50 Hz and plots the FFT spectra. This shows correct recovery above the Nyquist rate and aliasing below it.

The coefficient computation and reconstruction are written as reusable helper functions (`fs_coeffs`, `fs_reconstruct_Ck`, `stem_fs_spectrum`) rather than built-in toolbox shortcuts.

## Run

Open `fourier_series_analysis.m` in MATLAB (R2016b or later, with the Signal Processing Toolbox for `square`) and run it. Three figure windows are produced.

## Tech

MATLAB · Fourier analysis · FFT · sampling theory
