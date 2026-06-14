% Parameters
Fs = 1000;          % Sampling frequency (Hz)
T = 1;              % Duration (seconds)
t = 0:1/Fs:T-1/Fs;  % Time vector

% Generate pulse train signal
pulse_period = 0.2;          % Pulse every 0.2 seconds (5 Hz)
d = 0:pulse_period:0.8;      % Pulse delays: [0, 0.2, 0.4, 0.6, 0.8]
pulse_width = 0.05;          % Pulse width: 50 ms
x = pulstran(t, d, @rectpuls, pulse_width); % Rectangular pulse train

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
title('Spectrogram of Pulse Train Signal');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT
N = length(x);                % Signal length
X = fft(x);                   % Compute FFT
X_shifted = fftshift(X);      % Shift FFT to center at zero
freq = (-N/2:N/2-1)*(Fs/N);   % Frequency vector for shifted FFT

subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of Pulse Train Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Adjust layout for better spacing
sgtitle('Pulse Train Signal Analysis');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

saveas(gcf, 'imgs/signal_pulstran.png');
