% Parameters
fs = 10000;         % Sampling frequency (Hz)
t = 0:1/fs:2;       % Time vector from 0 to 2 s

% Generate VCO signal with sawtooth control
x = vco(sawtooth(2*pi*t, 0.75), [0.1 0.4]*fs, fs);

% Create single figure with subplots
figure;

% Subplot 1: Time-domain signal
subplot(3, 1, 1);
plot(t, x);
title('Time-Domain VCO Signal (Sawtooth Control)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Subplot 2: Spectrogram
subplot(3, 1, 2);
spectrogram(x, 256, 240, 256, fs, 'xaxis');
title('Spectrogram of VCO Signal (Sawtooth Control)');
colormap('jet');
colorbar;

% Subplot 3: Shifted FFT
N = length(x);
X = fft(x);
X_shifted = fftshift(X);
freq = (-N/2:N/2-1)*(fs/N);
subplot(3, 1, 3);
plot(freq, abs(X_shifted));
title('Shifted FFT of VCO Signal (Sawtooth Control)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
xlim([-5000, 5000]); % Zoom to capture relevant frequencies

% Adjust layout for better spacing
sgtitle('VCO Signal Analysis (Sawtooth Control)');
set(gcf, 'Position', [100, 100, 800, 800]); % Adjust figure size

% Save figure
saveas(gcf, 'imgs/signal_vco_sawtooth.png');