% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = 0:1/Fs:T-1/Fs;  % Time vector
f0 = 100;           % Fundamental frequency (Hz)

% Generate signals
x = square(2 * pi * f0 * t); % Square wave

% Figure 1: Square Wave Analysis
figure;

% Subplot 1: Time-domain signal (Square)
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Square Wave');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 100e-3]);
grid on;

% Subplot 2: Spectrogram (Square)
subplot(3, 1, 2);
spectrogram(x, 128, 120, 128, Fs, 'xaxis');
title('Spectrogram of Square Wave');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT (Square)
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(Fs/N);

subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Square Wave');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-400, 400]); % Zoom for clarity

% Adjust layout
sgtitle('Square Wave Analysis');
set(gcf, 'Position', [100, 100, 800, 800]);

saveas(gcf, 'imgs/signal_square.png');
