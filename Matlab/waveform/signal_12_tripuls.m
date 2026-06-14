% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = -T/2:1/Fs:T/2-1/Fs;  % Time vector centered at 0
f0 = 100;           % Reference frequency for scaling (Hz)

% Generate triangular pulse signal using tripuls function
w = 1/f0;           % Pulse width (0.01 s)
x = tripuls(t, w);  % Symmetric triangular pulse

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal (Triangular Pulse with stem)
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain Triangular Pulse Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([-20e-3, 20e-3]);
grid on;

% Subplot 2: Spectrogram (Triangular Pulse)
subplot(3, 1, 2);
spectrogram(x, 256, 240, 256, Fs, 'xaxis');
title('Spectrogram of Triangular Pulse Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT (Triangular Pulse)
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(Fs/N);
subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Triangular Pulse Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-500, 500]); % Zoom for clarity

% Adjust layout for better spacing
sgtitle('Triangular Pulse Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

% Save figure
saveas(gcf, 'imgs/signal_tripuls.png');