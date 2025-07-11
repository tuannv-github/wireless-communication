
pusch_channel = PuschChannel();
interference = Interference(25, -1e6);

tx_bits = randi([0 1], 50000, 1);
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

% Calculate and print throughput
throughput = length(tx_bits) / time_s;  % bits per second
fprintf('Throughput:\n');
fprintf('%.12f bits/second\n', throughput);
fprintf('%.12f kbits/second\n', throughput/1000);
fprintf('%.12f Mbits/second\n', throughput/1000000);
fprintf('%.12f Gbits/second\n', throughput/1000000000);
