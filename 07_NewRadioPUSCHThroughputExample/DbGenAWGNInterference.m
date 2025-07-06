clear;
clc;

mkdir('datasets/FccIQ/synthetic/ground_truth/singletone/rgb');
mkdir('datasets/FccIQ/synthetic/test/singletone/rgb');

SNRs = [5, 15];
SIRs = [5, 15];
MCSs = [14, 28];
FRQs = [-1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6];

% SNRs = [2, 4, 6, 8, 10, 14, 18, 22, 30, 38, 46];
% SIRs = [2, 4, 6, 8, 10, 14, 18, 22, 30, 38, 46];
% MCSs = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
% FRQ_min = -2e6;
% FRQ_max = 2e6;
% STEP = 1e6;
% FRQs = FRQ_min:STEP:FRQ_max;
% fprintf('FRQs: ');
% fprintf('%d ', FRQs);
% fprintf('\n');

% SIRs = [2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50];
% SNRs = [2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50];
% MCSs = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
% FRQ_min = -2e6;
% FRQ_max = 2e6;
% STEP = 1e6;
% FRQs = FRQ_min:STEP:FRQ_max;
% fprintf('FRQs: ');
% fprintf('%d ', FRQs);
% fprintf('\n');

function saveGrids(database, path, interference_type)
    % Define color limits for consistent visualization
    color_limits = [0, 2];
    
    for i = 1:length(database)
        item = database(i);
        fprintf('Saving SNR: %d, SIR: %d, MCS: %d, FRQ: %d\n', item.SNR, item.SIR, item.MCS, item.FRQ);
        % Extract NI data for each transmission
        for grid_idx = 1:length(item.grids)
            slot = item.grids(grid_idx).Slot;

            RGB_FILE_NAME_RX = sprintf('%s/test/%s/rgb/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_RxGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            PNG_FILE_NAME_RX = sprintf('%s/test/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_RxGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            MAT_FILE_NAME_RX = sprintf('%s/test/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_RxGrid.mat', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);

            RGB_FILE_NAME_NI = sprintf('%s/test/%s/rgb/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_NoiseInterferenceGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            PNG_FILE_NAME_NI = sprintf('%s/test/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_NoiseInterferenceGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            MAT_FILE_NAME_NI = sprintf('%s/test/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_NoiseInterferenceGrid.mat', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            
            RGB_FILE_NAME_INTERFERENCE = sprintf('%s/ground_truth/%s/rgb/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_InterferenceGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            PNG_FILE_NAME_INTERFERENCE = sprintf('%s/ground_truth/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_InterferenceGrid.png', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);
            MAT_FILE_NAME_INTERFERENCE = sprintf('%s/ground_truth/%s/SNR_%03d_SIR_%03d_MCS_%03d_FRQ_%09d_Slot_%03d_InterferenceGrid.mat', path, interference_type, item.SNR, item.SIR, item.MCS, item.FRQ, slot);

            % Check if files already exist
            if exist(PNG_FILE_NAME_RX, 'file') && exist(MAT_FILE_NAME_RX, 'file') && ...
               exist(PNG_FILE_NAME_NI, 'file') && exist(MAT_FILE_NAME_NI, 'file') && ...
               exist(PNG_FILE_NAME_INTERFERENCE, 'file') && exist(MAT_FILE_NAME_INTERFERENCE, 'file')
                fprintf('Files already exist for SNR: %d, SIR: %d, MCS: %d, FRQ: %d, slot: %d, skipping...\n', ...
                    item.SNR, item.SIR, item.MCS, item.FRQ, slot);
                continue;
            end

            rxGrid = item.grids(grid_idx).Rx;
            noiseInterferenceGrid = item.grids(grid_idx).NoiseInterference;
            interferenceGrid = item.grids(grid_idx).Interference;

            rxGrid_I = real(rxGrid);
            rxGrid_Q = imag(rxGrid);
            rxGrid_zero = zeros(size(rxGrid_I));
            rxGrid_IQ = cat(3, rxGrid_I, rxGrid_Q, rxGrid_zero);

            noiseInterferenceGrid_I = real(noiseInterferenceGrid);
            noiseInterferenceGrid_Q = imag(noiseInterferenceGrid);
            noiseInterferenceGrid_zero = zeros(size(noiseInterferenceGrid_I));
            noiseInterferenceGrid_IQ = cat(3, noiseInterferenceGrid_I, noiseInterferenceGrid_Q, noiseInterferenceGrid_zero);

            interferenceGrid_I = real(interferenceGrid);
            interferenceGrid_Q = imag(interferenceGrid);
            interferenceGrid_zero = zeros(size(interferenceGrid_I));
            interferenceGrid_IQ = cat(3, interferenceGrid_I, interferenceGrid_Q, interferenceGrid_zero);

            % Save IQ data as PNG and MAT files
            imwrite(rxGrid_IQ, PNG_FILE_NAME_RX);
            save(MAT_FILE_NAME_RX, 'rxGrid_IQ');

            imwrite(noiseInterferenceGrid_IQ, PNG_FILE_NAME_NI);
            save(MAT_FILE_NAME_NI, 'noiseInterferenceGrid_IQ');

            imwrite(interferenceGrid_IQ, PNG_FILE_NAME_INTERFERENCE);
            save(MAT_FILE_NAME_INTERFERENCE, 'interferenceGrid_IQ');

            % Create and save RGB visualization for RxGrid
            h = figure('Visible', 'off');
            imagesc(abs(rxGrid(:,:,1)));
            clim(color_limits);
            colorbar;
            title('Rx Grid Magnitude');
            xlabel('OFDM Symbols');
            ylabel('Subcarriers');
            saveas(h, RGB_FILE_NAME_RX);
            close(h);

            % Create and save RGB visualization for NoiseInterferenceGrid
            h = figure('Visible', 'off');
            imagesc(abs(noiseInterferenceGrid(:,:,1)));
            clim(color_limits);
            colorbar;
            title('Noise Interference Grid Magnitude');
            xlabel('OFDM Symbols');
            ylabel('Subcarriers');
            saveas(h, RGB_FILE_NAME_NI);
            close(h);

            % Create and save RGB visualization for InterferenceGrid
            h = figure('Visible', 'off');
            imagesc(abs(interferenceGrid(:,:,1)));
            clim(color_limits);
            colorbar;
            title('Interference Grid Magnitude');
            xlabel('OFDM Symbols');
            ylabel('Subcarriers');
            saveas(h, RGB_FILE_NAME_INTERFERENCE);
            close(h);
        end
    end
end

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
database = struct('SNR', {}, 'SIR', {}, 'MCS', {}, 'FRQ', {}, 'grids', {});

% Create parameter combinations for SIR, MCS, FRQ
[SIR_grid, MCS_grid, FRQ_grid] = ndgrid(SIRs, MCSs, FRQs);
params = [SIR_grid(:), MCS_grid(:), FRQ_grid(:)];

% Pre-allocate results array to track success
success_results = zeros(length(params), 1);

for SNR = SNRs
    database = [];
    parfor idx = 1:length(params)
        SIR = params(idx,1);
        MCS = params(idx,2);
        FRQ = params(idx,3);
        
        tx_bits = randi([0 1], 100000, 1);

        pusch_channel = PuschChannel();
        pusch_channel.simParameters.MCSIndex = MCS;
        pusch_channel.simParameters.SNR = SNR;
        pusch_channel.simParameters.Plot = false;
        pusch_channel.simParameters.DisplaySimulationInformation = false;
        pusch_channel.simParameters.Plot = false;

        interference = InterferenceSingletone(SIR);
        interference.toneFreq = FRQ;

        [rx_bits, time_s] = pusch_channel.tranceiver(tx_bits, interference);

        grids = pusch_channel.grids;
        if length(grids) > 0
            % Process noise interference data
            fprintf('Noise interference data available: %d samples\n', length(grids));
        end

        % Check if lengths match
        if length(tx_bits) ~= length(rx_bits)
            fprintf('Not all bits were received, received %d bits\n', length(rx_bits));
            success_results(idx) = 0;
            continue;
        end

        % Compare transmitted and received bits
        bit_errors = sum(tx_bits ~= rx_bits);
        if bit_errors ~= 0
            fprintf('Failed: %d bit errors detected\n', bit_errors);
            success_results(idx) = 0;
            continue;
        end

        % Mark as successful
        success_results(idx) = 1;
        fprintf('Success: All bits transmitted correctly\n');

        % Calculate throughput
        throughput = length(tx_bits) / time_s;  % bits per second
        fprintf('Throughput:\n');
        fprintf('%.12f bits/second\n', throughput);
        fprintf('%.12f kbits/second\n', throughput/1000);
        fprintf('%.12f Mbits/second\n', throughput/1000000);
        fprintf('%.12f Gbits/second\n', throughput/1000000000);
        fprintf('--------------------------------\n');

        row = struct('SNR', SNR, 'SIR', SIR, 'MCS', MCS, "FRQ", FRQ, "grids", grids);
        database = [database; row];
    end
    saveGrids(database, sprintf('datasets/FccIQ/synthetic/'), "singletone");
end
