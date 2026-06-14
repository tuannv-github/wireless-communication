% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = 0:1/Fs:T-1/Fs;  % Time vector
f0 = 100;           % Single-tone frequency (Hz)
A = 1;              % Amplitude

% Generate single-tone signal (sine wave)
x = A * sin(2 * pi * f0 * t);

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 100e-3]);
grid on;

% Subplot 2: Spectrogram
subplot(3, 1, 2);
spectrogram(x, 128, 120, 128, Fs, 'xaxis');
title('Spectrogram of Single-Tone Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT
N = length(x);                % Signal length
X = fft(x);                   % Compute FFT
X_shifted = fftshift(X);      % Shift FFT to center at zero
freq = (-N/2:N/2-1)*(Fs/N);   % Frequency vector for shifted FFT

subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Single-Tone Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Adjust layout for better spacing
sgtitle('Single-Tone Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

saveas(gcf, 'imgs/signal_sine.png');
