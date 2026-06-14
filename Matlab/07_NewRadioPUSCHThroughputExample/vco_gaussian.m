% VCO-controlled single tone with sawtooth control and spectrogram
clear all;
close all;

% Time vector
fs = 1e3; % Sampling frequency (Hz)
t = 10:1/fs:100; % Time range for signal

% VCO parameters
f0 = 0; % Center frequency (Hz)
Kv = 1; % VCO gain (Hz/V)

% Sawtooth control voltage
f_control = 1; % Frequency of sawtooth wave (Hz)
V_control = sin(2 * pi * f_control * t); % Sawtooth wave, amplitude [-1, 1]

% Calculate time-varying output frequency from VCO
f_vco = f0 + Kv * V_control;

% Generate single tone with time-varying frequency
sine_out = sin(2*pi*f_vco.*t); % Generate sine wave

% Plot time-domain signal
figure;
subplot(3,1,1);
plot(t, sine_out, 'b-', 'LineWidth', 2);
grid on;
title(['Single Tone with Sawtooth-Controlled VCO (f_0 = ', num2str(f0), ' Hz)']);
xlabel('Time (s)');
ylabel('Amplitude');

% Plot control voltage
subplot(3,1,2);
plot(t, V_control, 'r-', 'LineWidth', 2);
grid on;
title(['Sawtooth Control Voltage (f = ', num2str(f_control), ' Hz)']);
xlabel('Time (s)');
ylabel('Voltage (V)');

% Compute and plot spectrogram
subplot(3,1,3);
window = hamming(256); % Window size for spectrogram
noverlap = 128; % Overlap between windows
nfft = 512; % Number of FFT points
spectrogram(sine_out, window, noverlap, nfft, fs, 'yaxis');
title('Spectrogram of Single Tone');
colormap('jet');
colorbar;