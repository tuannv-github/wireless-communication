function [ret, row] = runSim(MCS, SNR, SIR, FRQ, tx_bits)
    ret = false;
    row = struct('MCS', MCS, 'SNR', SNR, 'SIR', SIR, 'FRQ', FRQ, 'Throughput', 0, 'RxGrids', []);

    fprintf('MCS: %d, SNR: %d, SIR: %d, FRQ: %d\n', MCS, SNR, SIR, FRQ);

    % Create local copies for parallel processing
    pusch_channel = PuschChannel();

    pusch_channel.simParameters.MCSIndex = MCS;
    pusch_channel.simParameters.SNR = SNR;
    interference = InterferenceSingletone(SIR);

    [rx_bits, time_s] = pusch_channel.tranceiver(tx_bits, interference);

    % Check if lengths match
    if length(tx_bits) ~= length(rx_bits)
        fprintf('Not all bits were received, received %d bits\n', length(rx_bits));
        return;
    end

    % Compare transmitted and received bits
    bit_errors = sum(tx_bits ~= rx_bits);
    if bit_errors ~= 0
        fprintf('Failed: %d bit errors detected\n', bit_errors);
        return;
    end

    fprintf('Success: All bits transmitted correctly\n');

    % Calculate throughput
    throughput = length(tx_bits) / time_s;  % bits per second
    fprintf('Throughput:\n');
    fprintf('%.12f bits/second\n', throughput);
    fprintf('%.12f kbits/second\n', throughput/1000);
    fprintf('%.12f Mbits/second\n', throughput/1000000);
    fprintf('%.12f Gbits/second\n', throughput/1000000000);
    fprintf('--------------------------------\n');

    row = struct('MCS', MCS, 'SNR', SNR, 'SIR', SIR, 'FRQ', FRQ, 'Throughput', throughput, 'RxGrids', pusch_channel.rxGrids);
    ret = true;
end
