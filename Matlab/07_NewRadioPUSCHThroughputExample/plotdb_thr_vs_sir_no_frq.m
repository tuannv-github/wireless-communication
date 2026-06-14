% Load the database
% load('database.mat');

% Extract unique values for each parameter
unique_MCSs = unique([database.MCS]);
unique_SNRs = unique([database.SNR]);
unique_SIRs = unique([database.SIR]);

% Define colors and markers
colors = {'b', 'r', 'g', 'm', 'c', 'k', [0.5 0.5 0], [0 0.5 0.5]};
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};

% Create all combinations of colors and markers
color_marker_combinations = {};
for i = 1:length(colors)
    for j = 1:length(markers)
        color_marker_combinations{end+1} = [colors{i}, markers{j}];
    end
end

% Create figure
figure('Name', 'Throughput vs SIR for Different Parameters');

% Plot for each MCS
for mcs_idx = 1:length(unique_MCSs)
    MCS = unique_MCSs(mcs_idx);
    
    % Create subplot for each MCS
    subplot(length(unique_MCSs), 1, mcs_idx);
    hold on;
    
    % Plot for each SNR
    for snr_idx = 1:length(unique_SNRs)
        SNR = unique_SNRs(snr_idx);
        
        % Initialize array to store averaged throughput values
        avg_throughput = zeros(length(unique_SIRs), 1);
        
        % For each SIR value
        for sir_idx = 1:length(unique_SIRs)
            SIR = unique_SIRs(sir_idx);
            
            % Get all data points for current MCS, SNR, and SIR
            mask = ([database.MCS] == MCS) & ([database.SNR] == SNR) & ([database.SIR] == SIR);
            throughput = [database.Throughput];
            current_throughput = throughput(mask);
            
            if ~isempty(current_throughput)
                % Average throughput over all FRQ values
                avg_throughput(sir_idx) = mean(current_throughput);
            end
        end
        
        % Plot averaged throughput vs SIR
        plot(unique_SIRs, avg_throughput/1e6, ...
            [color_marker_combinations{snr_idx}, '-'], ...
            'DisplayName', sprintf('SNR=%d', SNR));
    end
    
    % Add labels and title
    xlabel('SIR (dB)');
    ylabel('Throughput (Mbps)');
    title(sprintf('MCS = %d', MCS));
    grid on;
    legend('show', 'Location', 'best');
end

% Adjust figure layout
set(gcf, 'Position', [100, 100, 800, 1200]);
