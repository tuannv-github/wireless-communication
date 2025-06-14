
pusch_channel = PuschChannel();

tx_bits = randi([0 1], 10000000, 1);
[rx_bits, time_s] = pusch_channel.tranceiver(tx_bits);

% Compare transmitted and received bits
bit_errors = sum(tx_bits ~= rx_bits);
if bit_errors == 0
    fprintf('Success: All bits transmitted correctly\n');
else
    fprintf('Failed: %d bit errors detected\n', bit_errors);
end

% Calculate and print throughput
throughput = length(tx_bits) / time_s;  % bits per second
fprintf('Throughput:\n');
fprintf('%.12f bits/second\n', throughput);
fprintf('%.12f kbits/second\n', throughput/1000);
fprintf('%.12f Mbits/second\n', throughput/1000000);
fprintf('%.12f Gbits/second\n', throughput/1000000000);
