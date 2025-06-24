clear;
clc;

database = [];

if ~exist('database', 'dir')
    mkdir('database');
end

% if exist('database', 'dir')
%     delete('database/*');
% end

MCSs = [15];
SNRs = [40];
SIRs = [15];
FRQs = [1e6];

% MCSs = [14, 28];
% SNRs = [0, 20, 40];
% SIRs = [0, 20, 40];
% FRQs = [-1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6];

% MCSs = [14, 28];
% SNRs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
% SIRs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
% FRQs = [-2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6];

% MCSs = [14, 28];
% SNRs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
% SIRs = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50];
% FRQs = [-2e6, -1.5e6, -1e6, -0.5e6, 0, 0.5e6, 1e6, 1.5e6, 2e6];

% MCSs = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];
% SNRs = [0, 2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50];
% SIRs = [0, 2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50];
% FRQ_min = -2e6;
% FRQ_max = 2e6;
% STEP = 0.1e6;
% FRQs = FRQ_min:STEP:FRQ_max;

tx_bits = randi([0 1], 100000, 1);

RUN_IN_PARALLEL = false;

if ~RUN_IN_PARALLEL
    for MCS = MCSs
        for SNR = SNRs
            for SIR = SIRs
                for FRQ = FRQs
                    row = runSim(MCS, SNR, SIR, FRQ, tx_bits);
                    database = [database; row];
                end
            end
        end
    end
else
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
    num_combinations = length(MCSs) * length(SNRs) * length(SIRs) * length(FRQs);
    database = struct('MCS', {}, 'SNR', {}, 'SIR', {}, 'FRQ', {}, 'Throughput', {}, 'RxGrids', {});

    % % Create parameter combinations
    % [MCS_grid, SNR_grid, SIR_grid, FRQ_grid] = ndgrid(MCSs, SNRs, SIRs, FRQs);
    % params = [MCS_grid(:), SNR_grid(:), SIR_grid(:), FRQ_grid(:)];

    for MCS = MCSs
        for SNR = SNRs
            for SIR = SIRs
                filename = sprintf('database/pusch_simulation_database_MCS_%d_SNR_%d_SIR_%d.mat', MCS, SNR, SIR);
                if exist(filename, 'file')
                    fprintf('Skipping MCS=%d, SNR=%d, SIR=%d - file already exists\n', MCS, SNR, SIR);
                    continue;
                end
                database = [];
                parfor i = 1:length(FRQs)
                    FRQ = FRQs(i);
                    [ret, row] = runSim(MCS, SNR, SIR, FRQ, tx_bits);
                    if ret
                        database = [database; row];
                    end
                end
                save(filename, 'database');
            end
        end
    end
end
