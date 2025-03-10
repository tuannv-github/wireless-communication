% MATLAB script to extract ToF from CSI using MUSIC algorithm

% System Parameters
fs = 20e6;           % Sampling frequency (20 MHz, typical for Wi-Fi)
delta_f = 312.5e3;   % Subcarrier spacing (e.g., 312.5 kHz for 802.11n)
Nsc = 64;            % Number of subcarriers
c = 3e8;             % Speed of light (m/s)
Ts = 1/fs;           % Sampling period

% Generate synthetic CSI with multipath
true_tofs = [50e-9, 100e-9, 150e-9];  % True ToFs in seconds (e.g., 50ns, 100ns, 150ns)
num_paths = length(true_tofs);
amplitudes = [1, 0.8, 0.6];          % Path amplitudes
noise_power = 0.01;                  % Noise power

% Frequency vector for subcarriers
f = (-(Nsc/2 - 1):Nsc/2) * delta_f;  % Subcarrier frequencies (centered)

% Generate CIR in frequency domain (CSI)
csi = zeros(1, Nsc);
for p = 1:num_paths
    csi = csi + amplitudes(p) * exp(-1j * 2 * pi * f * true_tofs(p));
end
csi = csi + sqrt(noise_power/2) * (randn(1, Nsc) + 1j * randn(1, Nsc));  % Add noise

% MUSIC Algorithm
% Step 1: Compute covariance matrix
R = (csi' * csi) / Nsc;  % Simplified for single snapshot; ideally use multiple snapshots

% Step 2: Eigendecomposition
[E, D] = eig(R);
eigenvalues = diag(D);
[eigenvalues, idx] = sort(eigenvalues, 'descend');  % Sort eigenvalues
E = E(:, idx);  % Sort eigenvectors

% Step 3: Estimate number of paths (e.g., using MDL or manual)
num_estimated_paths = num_paths;  % Assume known for simplicity

% Step 4: Noise subspace ( eigenvectors corresponding to smallest eigenvalues)
En = E(:, num_estimated_paths+1:end);

% Step 5: Define delay search range and MUSIC spectrum
tau_max = 200e-9;  % Max ToF to search (200 ns)
tau_step = 1e-9;   % Resolution (1 ns)
taus = 0:tau_step:tau_max;  % Delay search vector
P_music = zeros(size(taus));

% Steering vector function
steering = @(tau) exp(-1j * 2 * pi * f' * tau);

% Compute MUSIC spectrum
for i = 1:length(taus)
    a = steering(taus(i));  % Steering vector for delay tau
    P_music(i) = 1 / (a' * (En * En') * a);  % MUSIC pseudo-spectrum
end

% Normalize spectrum
P_music = abs(P_music) / max(abs(P_music));

% Step 6: Peak detection to find ToFs
[~, peaks] = findpeaks(P_music, 'MinPeakHeight', 0.1, 'SortStr', 'descend', 'NPeaks', num_paths);
estimated_tofs = taus(peaks);

% Convert to distance (optional)
estimated_distances = estimated_tofs * c;

% Display results
disp('True ToFs (ns):');
disp(true_tofs * 1e9);
disp('Estimated ToFs (ns):');
disp(estimated_tofs * 1e9);

% Plot MUSIC spectrum
figure;
plot(taus * 1e9, P_music, 'b-', 'LineWidth', 1.5);
hold on;
stem(true_tofs * 1e9, ones(size(true_tofs)), 'r--', 'LineWidth', 1.5, 'Marker', 'o');
stem(estimated_tofs * 1e9, ones(size(estimated_tofs)), 'g--', 'LineWidth', 1.5, 'Marker', 'x');
xlabel('Delay (ns)');
ylabel('MUSIC Spectrum');
title('MUSIC Algorithm for ToF Estimation');
legend('MUSIC Spectrum', 'True ToFs', 'Estimated ToFs');
grid on;

% Plot CSI (magnitude)
figure;
plot(f/1e6, abs(csi), 'b-', 'LineWidth', 1.5);
xlabel('Frequency (MHz)');
ylabel('CSI Magnitude');
title('Synthetic CSI Magnitude');
grid on;