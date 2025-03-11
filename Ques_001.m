clc;
clear;
close all;

% Parameters
M_values = [2, 4, 8]; % Modulation orders (BPSK, QPSK, 8-PSK)
SNR_range = 0:15; % SNR range in dB
numBits = 100000; % Number of bits for simulation

% Initialize BER results
BER = zeros(length(M_values), length(SNR_range));

% Loop over modulation orders
for m_idx = 1:length(M_values)
    M = M_values(m_idx);
    bps = log2(M); % Bits per symbol
    
    % Generate random binary data
    data = randi([0 1], 1, numBits);
    
    % Padding to ensure data length is a multiple of bps
    remainder = mod(numBits, bps);
    if remainder ~= 0
        paddingBits = zeros(1, bps - remainder);
        data = [data paddingBits];
    end
    
    % Convert bits to symbols
    reshapedData = reshape(data, [], bps);
    symbols = bi2de(reshapedData, 'left-msb');
    
    % PSK modulation
    modulatedSignal = pskmod(symbols, M, 0);
    
    % Loop over SNR values
    for snr_idx = 1:length(SNR_range)
        snr = SNR_range(snr_idx);
        
        % Add AWGN noise
        noisySignal = awgn(modulatedSignal, snr, 'measured');
        
        % PSK demodulation
        demodulatedSymbols = pskdemod(noisySignal, M, 0);
        
        % Convert symbols back to bits
        demodulatedBits = de2bi(demodulatedSymbols, bps, 'left-msb');
        receivedBits = reshape(demodulatedBits.', 1, []);
        
        % Remove padding
        receivedBits = receivedBits(1:numBits);
        
        % Calculate BER
        [~, ber] = biterr(data(1:numBits), receivedBits);
        BER(m_idx, snr_idx) = ber;
    end
end

% Plot BER vs SNR
figure;
semilogy(SNR_range, BER(1, :), 'b-o', 'LineWidth', 2);
hold on;
semilogy(SNR_range, BER(2, :), 'r-s', 'LineWidth', 2);
semilogy(SNR_range, BER(3, :), 'g-d', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for PSK Modulation');
legend('BPSK (M=2)', 'QPSK (M=4)', '8-PSK (M=8)');
grid on;