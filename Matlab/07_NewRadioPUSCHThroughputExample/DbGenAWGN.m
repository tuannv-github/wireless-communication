clear;
clc;

% if exist('datasets/FccIQ/synthetic/train/good', 'dir')
%     fprintf('Deleting existing directory\n');
%     rmdir('datasets/FccIQ/synthetic/train/good', 's');
% end
fprintf('Creating new directory\n');
mkdir('datasets/FccIQ/synthetic/train/good/rgb');

% SNRs = [10, 20, 40];
% MCSs = [1, 14, 28];

SNRs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
MCSs = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];

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

    % Verify poonrChannelEstimatel size
    pool = gcp;
    fprintf('Parallel pool started with %d workers\n', pool.NumWorkers);
end

% Pre-allocate database array
num_combinations = length(MCSs) * length(SNRs);
database = struct('SNR', {}, 'MCS', {}, 'grids', {});

% Create parameter combinations
[SNR, MCS] = ndgrid(SNRs, MCSs);
params = [SNR(:), MCS(:)];

% Pre-allocate results array to track success
success_results = zeros(length(params), 1);

parfor i = 1:length(params)
    SNR = params(i,1);
    MCS = params(i,2);
    crc_count = 10000;

    pusch_channel = PuschChannel();
    pusch_channel.simParameters.MCSIndex = MCS;
    pusch_channel.simParameters.SNR = SNR;
    pusch_channel.simParameters.Plot = false;
    pusch_channel.simParameters.DisplaySimulationInformation = false;
    pusch_channel.simParameters.Plot = false;

    ret = pusch_channel.tranceiver(crc_count);
    success_results(i) = ret;
    if ~ret
        fprintf('Failed: %d\n', i);
        continue;
    end

    grids = pusch_channel.grids;
    if length(grids) > 0
        % Process noise interference data
        fprintf('Noise interference data available: %d samples\n', length(grids));
    end

    row = struct('SNR', SNR, 'MCS', MCS, 'grids', pusch_channel.grids);
    database = [database; row];
end

% Sort database by SNR and MCS
[~, sort_idx] = sortrows([cell2mat({database.SNR})', cell2mat({database.MCS})'], [1, 2]);
database = database(sort_idx);
fprintf('Database sorted by SNR and MCS\n');

% Sort success_results by SNR and MCS
[~, sort_indices] = sortrows([params(:,1), params(:,2)], [1, 2]); % Sort by SNR first, then MCS
success_results = success_results(sort_indices);

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
        SNR = params(idx,1);
        MCS = params(idx,2);
        fprintf('SNR: %.1f dB, MCS: %d\n', SNR, MCS);
    end
    fprintf('--------------------------------\n');
    fprintf('Total failed cases: %d\n', length(failed_indices));
else
    fprintf('\nAll transmissions were successful!\n');
end

% Save successful database entries to files
fprintf('\nSaving successful database entries...\n');

% Pre-define base paths to avoid repeated string concatenation
base_path = 'datasets/FccIQ/synthetic/train/good';
rgb_path = fullfile(base_path, 'rgb');

% Save the entire database to a single .mat file
database_file_path = fullfile(base_path, 'database.mat');
fprintf('Saving complete database to: %s\n', database_file_path);
save(database_file_path, 'database', '-v7.3');
fprintf('Database saved successfully!\n');

% % Pre-define color limits for consistency
% color_limits = [0, 2];

% for i = 1:length(database)
%     item = database(i);
%     fprintf('Saving SNR: %d, MCS: %d\n', item.SNR, item.MCS);
    
%     % Extract NI data for each transmission
%     for grid_idx = 1:length(item.grids)
%         slot = item.grids(grid_idx).Slot;
        
%         % Generate file names once
%         file_suffix = sprintf('SNR_%d_MCS_%d_Slot_%d_GridIdx_%d', item.SNR, item.MCS, slot, grid_idx);
        
%         PNG_FILE_NAME_RX = fullfile(base_path, [file_suffix '_RxGrid.png']);
%         MAT_FILE_NAME_RX = fullfile(base_path, [file_suffix '_RxGrid.mat']);
%         RGB_FILE_NAME_RX = fullfile(rgb_path, [file_suffix '_RxGrid.png']);
        
%         PNG_FILE_NAME_NI = fullfile(base_path, [file_suffix '_NoiseInterferenceGrid.png']);
%         MAT_FILE_NAME_NI = fullfile(base_path, [file_suffix '_NoiseInterferenceGrid.mat']);
%         RGB_FILE_NAME_NI = fullfile(rgb_path, [file_suffix '_NoiseInterferenceGrid.png']);
        
%         % Check if all files already exist
%         all_files_exist = all([exist(PNG_FILE_NAME_RX, 'file'), exist(MAT_FILE_NAME_RX, 'file'), ...
%                               exist(RGB_FILE_NAME_RX, 'file'), exist(PNG_FILE_NAME_NI, 'file'), ...
%                               exist(MAT_FILE_NAME_NI, 'file'), exist(RGB_FILE_NAME_NI, 'file')]);
        
%         if all_files_exist
%             fprintf('Files already exist for SNR: %d, MCS: %d, slot: %d, grid_idx: %d, skipping...\n', ...
%                 item.SNR, item.MCS, slot, grid_idx);
%             continue;
%         end
        
%         % Extract grid data
%         rxGrid = item.grids(grid_idx).Rx;
%         noiseInterferenceGrid = item.grids(grid_idx).NoiseInterference;
        
%         % Process RxGrid
%         rxGrid_I = real(rxGrid);
%         rxGrid_Q = imag(rxGrid);
%         rxGrid_IQ = cat(3, rxGrid_I, rxGrid_Q, zeros(size(rxGrid_I)));
        
%         imwrite(rxGrid_IQ, PNG_FILE_NAME_RX);
%         save(MAT_FILE_NAME_RX, 'rxGrid_IQ');

%         % Create and save RGB visualization for RxGrid
%         h = figure('Visible', 'off');
%         imagesc(abs(rxGrid(:,:,1)));
%         clim(color_limits);
%         colorbar;
%         title('Rx Grid Magnitude');
%         xlabel('OFDM Symbols');
%         ylabel('Subcarriers');
%         saveas(h, RGB_FILE_NAME_RX);
%         close(h);

%         % Process NoiseInterferenceGrid
%         noiseInterferenceGrid_I = real(noiseInterferenceGrid);
%         noiseInterferenceGrid_Q = imag(noiseInterferenceGrid);
%         noiseInterferenceGrid_IQ = cat(3, noiseInterferenceGrid_I, noiseInterferenceGrid_Q, zeros(size(noiseInterferenceGrid_I)));
        
%         imwrite(noiseInterferenceGrid_IQ, PNG_FILE_NAME_NI);
%         save(MAT_FILE_NAME_NI, 'noiseInterferenceGrid_IQ');
        
%         % Create and save RGB visualization for NoiseInterferenceGrid
%         h = figure('Visible', 'off');
%         imagesc(abs(noiseInterferenceGrid(:,:,1)));
%         clim(color_limits);
%         colorbar;
%         title('Noise Interference Grid Magnitude');
%         xlabel('OFDM Symbols');
%         ylabel('Subcarriers');
%         saveas(h, RGB_FILE_NAME_NI);
%         close(h);
%     end
% end
