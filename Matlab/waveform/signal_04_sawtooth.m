% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = 0:1/Fs:T-1/Fs;  % Time vector
f0 = 100;           % Fundamental frequency (Hz)

% Generate signals
x = sawtooth(2 * pi * f0 * t); % Sawtooth wave

% Figure 2: Sawtooth Wave Analysis
figure;

% Subplot 1: Time-domain signal (Sawtooth)
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Sawtooth Wave');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 100e-3]);
grid on;

% Subplot 2: Spectrogram (Sawtooth)
subplot(3, 1, 2);
spectrogram(x, 128, 120, 128, Fs, 'xaxis');
title('Spectrogram of Sawtooth Wave');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT (Sawtooth)
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(Fs/N);
subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Sawtooth Wave');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-400, 400]); % Zoom for clarity

% Adjust layout
sgtitle('Sawtooth Wave Analysis');
set(gcf, 'Position', [900, 100, 800, 800]); % Offset second figure

saveas(gcf, 'imgs/signal_sawtooth.png');
