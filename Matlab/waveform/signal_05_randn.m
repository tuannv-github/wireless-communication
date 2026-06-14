% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = 0:1/Fs:T-1/Fs;  % Time vector

% Generate random noise signal
x = randn(1, length(t)); % Gaussian white noise

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Subplot 2: Spectrogram
subplot(3, 1, 2);
spectrogram(x, 128, 120, 128, Fs, 'xaxis');
title('Spectrogram of Random Noise Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT
N = length(x);                % Signal length
X = fft(x);                   % Compute FFT
X_shifted = fftshift(X);      % Shift FFT to center at zero
freq = (-N/2:N/2-1)*(Fs/N);   % Frequency vector for shifted FFT

subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Random Noise Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Adjust layout for better spacing
sgtitle('Random Noise Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

saveas(gcf, 'imgs/signal_randn.png');
