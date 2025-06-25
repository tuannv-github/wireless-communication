clear;
clc;

if ~exist('dataset/FccIQ/good', 'dir')
    mkdir('dataset/FccIQ/good');
end

if exist('dataset/FccIQ/good', 'dir')
    delete('dataset/FccIQ/good/*');
end

if ~exist('dataset/FccIQ/good/rgb', 'dir')
    mkdir('dataset/FccIQ/good/rgb');
end

if exist('dataset/FccIQ/good/rgb', 'dir')
    delete('dataset/FccIQ/good/rgb/*');
end

% MCSs = [15];
% SNRs = [40];

% MCSs = [14, 28];
% SNRs = [0, 20, 40];

% MCSs = [14, 28];
% SNRs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50];

% MCSs = [14, 28];
% SNRs = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50];

MCSs = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
SNRs = [2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50];

fprintf('Running in parallel\n');

% Create a parallel pool if it doesn't exist
if isempty(gcp('nocreate'))
    % Get the cluster profile
    c = parcluster('Processes');

    numCores = feature('numcores');
    fprintf('Logical cores: %d\n', numCores);

    % Set NumWorkers to the desired number
    c.NumWorkers = numCores; % Adjust based on your hardware/license

    % Save the updated profile
    saveProfile(c);

    % Start parpool with the desired number of workers
    parpool('Processes', numCores);

    % Verify pool size
    pool = gcp;
    fprintf('Parallel pool started with %d workers\n', pool.NumWorkers);
end

% Pre-allocate database array
num_combinations = length(MCSs) * length(SNRs);
database = struct('MCS', {}, 'SNR', {}, 'NIs', {});

% Create parameter combinations
[MCS_grid, SNR_grid] = ndgrid(MCSs, SNRs);
params = [MCS_grid(:), SNR_grid(:)];

% Pre-allocate results array to track success
success_results = zeros(length(params), 1);

parfor i = 1:length(params)
    MCS = params(i,1);
    SNR = params(i,2);
    tx_bits = randi([0 1], 100000, 1);

    pusch_channel = PuschChannel();
    pusch_channel.simParameters.MCSIndex = MCS;
    pusch_channel.simParameters.SNR = SNR;
    pusch_channel.simParameters.Plot = false;
    pusch_channel.simParameters.DisplaySimulationInformation = false;
    pusch_channel.simParameters.Plot = false;

    [rx_bits, time_s] = pusch_channel.tranceiver(tx_bits);

    NIs = pusch_channel.NIs;
    if length(NIs) > 0
        % Process noise interference data
        fprintf('Noise interference data available: %d samples\n', length(NIs));
    end

    % Check if lengths match
    if length(tx_bits) ~= length(rx_bits)
        fprintf('Not all bits were received, received %d bits\n', length(rx_bits));
        success_results(i) = 0;
        continue;
    end

    % Compare transmitted and received bits
    bit_errors = sum(tx_bits ~= rx_bits);
    if bit_errors ~= 0
        fprintf('Failed: %d bit errors detected\n', bit_errors);
        success_results(i) = 0;
        continue;
    end

    % Mark as successful
    success_results(i) = 1;
    fprintf('Success: All bits transmitted correctly\n');

    % Calculate throughput
    throughput = length(tx_bits) / time_s;  % bits per second
    fprintf('Throughput:\n');
    fprintf('%.12f bits/second\n', throughput);
    fprintf('%.12f kbits/second\n', throughput/1000);
    fprintf('%.12f Mbits/second\n', throughput/1000000);
    fprintf('%.12f Gbits/second\n', throughput/1000000000);
    fprintf('--------------------------------\n');

    row = struct('MCS', MCS, 'SNR', SNR, 'NIs', pusch_channel.NIs);
    database = [database; row];
end

% Calculate total successful transmissions after parfor loop
successful_transmissions = sum(success_results);
fprintf('Total successful transmissions: %d out of %d attempts\n', successful_transmissions, length(params));

% Print failed cases
failed_indices = find(success_results == 0);
if ~isempty(failed_indices)
    fprintf('\nFailed cases:\n');
    fprintf('--------------------------------\n');
    for i = 1:length(failed_indices)
        idx = failed_indices(i);
        MCS = params(idx,1);
        SNR = params(idx,2);
        fprintf('MCS: %d, SNR: %.1f dB\n', MCS, SNR);
    end
    fprintf('--------------------------------\n');
    fprintf('Total failed cases: %d\n', length(failed_indices));
else
    fprintf('\nAll transmissions were successful!\n');
end

% Save successful database entries to files
fprintf('\nSaving successful database entries...\n');
successful_count = 0;

for i = 1:length(database)
    item = database(i);
    fprintf('Saving MCS: %d, SNR: %d\n', item.MCS, item.SNR);
    % Extract NI data for each transmission
    for ni_idx = 1:length(item.NIs)
        slot = item.NIs(ni_idx).Slot;
        NI = item.NIs(ni_idx).NI;
        
        NI_I = real(NI);
        NI_Q = imag(NI);
        NI_zero = zeros(size(NI_I));
        NI_IQ = cat(3, NI_I, NI_Q, NI_zero);
        % Save H_IQ matrix as image without plotting

        imwrite(mat2gray(NI_IQ(:,:,1)), sprintf('dataset/FccIQ/good/MCS_%d_SNR_%d_NI_%d_slot_%d.png', item.MCS, item.SNR, ni_idx, slot));
        save(sprintf('dataset/FccIQ/good/MCS_%d_SNR_%d_NI_%d_slot_%d.mat', item.MCS, item.SNR, ni_idx, slot), 'NI_IQ');

        h = figure('Visible', 'off');
        imagesc(abs(NI(:,:,1)));
        clim([0, 2]);
        colorbar;
        title('Noise Interference Grid Magnitude');
        xlabel('OFDM Symbols');
        ylabel('Subcarriers');
        saveas(h, sprintf('dataset/FccIQ/good/rgb/MCS_%d_SNR_%d_NI_%d_slot_%d.png', item.MCS, item.SNR, ni_idx, slot));
        close(h);
    end
end

fprintf('Successfully saved %d database entries to dataset/FccIQ/good/\n', successful_count);
