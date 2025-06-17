% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = -T/2:1/Fs:T/2-1/Fs;  % Time vector centered at 0
f0 = 100;           % Reference frequency for scaling (Hz)

% Generate sinc signal using MATLAB's sinc function
x = sinc(2 * f0 * t); % sinc(πx) where x = 2*f0*t, so πx = 2πf0t

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal (Sinc)
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Sinc Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Subplot 2: Spectrogram (Sinc)
subplot(3, 1, 2);
spectrogram(x, 256, 240, 256, Fs, 'xaxis');
title('Spectrogram of Sinc Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT (Sinc)
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(Fs/N);
subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Sinc Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-500, 500]); % Zoom for clarity

% Adjust layout for better spacing
sgtitle('Sinc Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

saveas(gcf, 'imgs/signal_sinc.png');
