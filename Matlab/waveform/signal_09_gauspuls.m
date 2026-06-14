% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = -T/2:1/Fs:T/2-1/Fs;  % Time vector centered at 0
f0 = 100;           % Reference frequency for scaling (Hz)

% Generate Gaussian pulse signal using gauspuls function
x = gauspuls(t, f0, 1); % Gaussian pulse with center frequency f0 and bandwidth 1

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal (Gaussian Pulse)
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Gaussian Pulse Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Subplot 2: Spectrogram (Gaussian Pulse)
subplot(3, 1, 2);
spectrogram(x, 256, 240, 256, Fs, 'xaxis');
title('Spectrogram of Gaussian Pulse Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT (Gaussian Pulse)
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(Fs/N);
subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Gaussian Pulse Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-500, 500]); % Zoom for clarity

% Adjust layout for better spacing
sgtitle('Gaussian Pulse Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

% Save figure
saveas(gcf, 'imgs/signal_gauspuls.png');