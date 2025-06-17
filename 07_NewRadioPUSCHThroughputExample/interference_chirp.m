% Parameters
Fs = 10e6; % Sampling frequency (10 MHz)
T = 0.5; % Duration (1 second)
t = 0:1/Fs:T-1/Fs; % Time vector (Δt = 1e-7 s)
f0 = 0; % Initial frequency (0 Hz)
f1 = 1e3; % Final frequency (1 MHz)
A = 1; % Amplitude

% Generate chirp signal (linear frequency sweep)
x = A * chirp(t, f0, T, f1, 'linear');

% Compute FFT
N = length(x); % Number of samples
X = fft(x)/N; % Normalized FFT
freq = (-N/2:N/2-1)*(Fs/N); % Frequency vector (centered at 0 Hz)
X_shifted = fftshift(X); % Shift FFT for negative frequencies
X_mag = abs(X_shifted); % Magnitude spectrum

% Apply 1 MHz bandwidth filter (low-pass, -1 MHz to 1 MHz)
f_cutoff = 0.5e6; % 1 MHz cutoff frequency
filter_mask = abs(freq) <= f_cutoff; % Keep frequencies |f| ≤ 1 MHz
X_filtered = X_shifted .* filter_mask; % Apply filter
X_filtered_unshifted = ifftshift(X_filtered); % Unshift for IFFT

% Inverse FFT to get filtered time-domain signal
x_filtered = ifft(X_filtered_unshifted)*N; % Inverse FFT, rescaled

% Plot results
figure;

% Subplot 1: Original chirp signal
subplot(4,1,1);
plot(t, x);
title('Chirp Signal (0 Hz to 1 MHz)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 50e-3]); % Zoom to 1 ms to show frequency increase
grid on;

% Subplot 2: Original and filtered frequency spectrum (with negative frequencies)
subplot(4,1,2);
plot(freq, X_mag, 'b', 'DisplayName', 'Original');
hold on;
title('Frequency Domain: Original vs Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10e3 10e3]);
legend;
grid on;

% Subplot 2: Original and filtered frequency spectrum (with negative frequencies)
subplot(4,1,3);
hold on;
plot(freq, abs(X_filtered), 'r', 'DisplayName', 'Filtered (1 MHz BW)');
title('Frequency Domain: Original vs Filtered Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10e3 10e3]);
legend;
grid on;

% Subplot 3: Original and filtered time-domain signal
subplot(4,1,4);
plot(t, x, 'b', 'DisplayName', 'Original');
hold on;
plot(t, real(x_filtered), 'r', 'DisplayName', 'Filtered');
title('Time Domain: Original vs Filtered Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 1]); % Zoom to 1 ms to show frequency increase
legend;
grid on;