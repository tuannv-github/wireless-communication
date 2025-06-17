% Parameters
fs = 44100;        % Sampling frequency (Hz)
T = 1;             % Duration (seconds)
t = 0:1/fs:T;      % Time vector
f0 = 100;          % Start frequency (Hz)
f1 = 1000;         % End frequency (Hz)

% Generate chirp signal
signal = chirp(t, f0, T, f1, 'linear');

% Compute FFT
N = length(signal);
X = fft(signal);
X_shifted = fftshift(X);
f = (-N/2:N/2-1)*(fs/N);

% Plotting
figure('Position', [100, 100, 800, 800]);

% Time-domain signal
subplot(3, 1, 1);
plot(t, signal);
title('Linear Chirp Signal (100 Hz to 1 kHz)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Spectrogram
subplot(3, 1, 2);
spectrogram(signal, 256, 128, 256, fs, 'yaxis');
title('Spectrogram');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colormap('jet');
colorbar;
% ylim([0 1500]); % Limit y-axis to focus on 0-1.5 kHz range

% FFT magnitude spectrum
subplot(3, 1, 3);
plot(f, abs(X_shifted)/N);
title('Magnitude Spectrum (FFT with fftshift)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-2000 2000]);

% Adjust layout
set(gcf, 'Color', 'white');