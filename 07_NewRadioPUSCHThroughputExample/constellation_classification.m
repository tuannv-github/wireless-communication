% PUSCH Modulation Scheme Detection from IQ Resource Grid
clear; clc;

% Simulate IQ resource grid (replace with actual data)
numSubcarriers = 12 * 4; % 4 RBs, 12 subcarriers per RB
numSymbols = 14; % 1 slot with 14 OFDM symbols
modOrder = 16; % Simulate 16-QAM (4: QPSK, 16: 16-QAM, 64: 64-QAM, 256: 256-QAM)
snr_dB = 20; % Signal-to-noise ratio in dB

% Generate random PUSCH data
bits = randi([0, 1], log2(modOrder) * numSubcarriers * (numSymbols-2), 1);
modSymbols = qammod(bits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);

% Create resource grid (complex IQ symbols)
iqGrid = complex(zeros(numSubcarriers, numSymbols));
puschIdx = 1;
for sym = 1:numSymbols
    if sym == 3 || sym == 10 % Assume DMRS in symbols 3 and 10
        iqGrid(:, sym) = randn(numSubcarriers, 1) + 1j * randn(numSubcarriers, 1); % DMRS
    else
        iqGrid(:, sym) = modSymbols(puschIdx:puschIdx+numSubcarriers-1);
        puschIdx = puschIdx + numSubcarriers;
    end
end

% Add AWGN noise
iqGrid = awgn(iqGrid, snr_dB, 'measured');

% Step 1: Extract PUSCH and DMRS symbols
puschSymbols = [];
for sym = 1:numSymbols
    if sym ~= 3 && sym ~= 10 % Exclude DMRS symbols
        puschSymbols = [puschSymbols; iqGrid(:, sym)];
    end
end

% Step 2: Channel estimation using DMRS (simplified)
dmrsSymbols = [iqGrid(:, 3); iqGrid(:, 10)]; % Received DMRS
idealDmrs = ones(size(dmrsSymbols)); % Assume ideal DMRS = 1 (simplified)
channelEst = mean(dmrsSymbols ./ idealDmrs); % Average channel estimate

% Step 3: Equalize PUSCH symbols
equalizedSymbols = puschSymbols / channelEst;

% Step 4: Plot constellation
figure;
plot(real(puschSymbols), imag(puschSymbols), 'o', 'MarkerSize', 4);
grid on;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');

% Step 5: Plot constellation
figure;
plot(real(equalizedSymbols), imag(equalizedSymbols), 'o', 'MarkerSize', 4);
grid on;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');