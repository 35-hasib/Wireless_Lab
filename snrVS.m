clc;
clear all;
close all;

%% Parameters
M = 8; % Modulation order (8-PSK)
bps = log2(M); % Bits per symbol

%% Input Text to Binary Conversion
txt1 = 'Information and communication engineering';
symbols = double(txt1); % Convert text to ASCII values
symbolToBitMapping = de2bi(symbols, 8, 'left-msb'); % Convert ASCII to binary

totNoBits = numel(symbolToBitMapping); % Total number of bits
inputReshapedBits = reshape(symbolToBitMapping, 1, totNoBits); % Reshape to 1D array

%% Padding
remainder = rem(totNoBits, bps);
if remainder == 0
    userPaddedData = inputReshapedBits;
else
    paddingBits = zeros(1, bps - remainder); % Add padding bits
    userPaddedData = [inputReshapedBits, paddingBits];
end

%% Modulation
reshapedUserPaddedData = reshape(userPaddedData, numel(userPaddedData)/bps, bps);
bitToSymbolMapping = bi2de(reshapedUserPaddedData, 'left-msb'); % Convert bits to symbols
modulatedSymbol = pskmod(bitToSymbolMapping, M); % 8-PSK modulation

%% Channel Simulation (SNR Variation)
SNR = [];
BER = [];

for snr = 0:15
    SNR = [SNR, snr]; % Store SNR values
    noisySymbols = awgn(modulatedSymbol, snr, 'measured'); % Add AWGN noise
    demodulatedSymbol = pskdemod(noisySymbols, M); % 8-PSK demodulation

    % Convert demodulated symbols back to bits
    demodulatedSymbolToBitMapping = de2bi(demodulatedSymbol, bps, 'left-msb');
    reshapedDemodulatedBits = reshape(demodulatedSymbolToBitMapping, 1, numel(demodulatedSymbolToBitMapping));

    % Remove padding bits
    demodulatedBitsWithoutPadding = reshapedDemodulatedBits(1:totNoBits);

    % Calculate Bit Error Rate (BER)
    [noe, ber] = biterr(inputReshapedBits, demodulatedBitsWithoutPadding);
    BER = [BER, ber]; % Store BER values

    % Reconstruct original text from demodulated bits
    txtBits = reshape(demodulatedBitsWithoutPadding, numel(demodulatedBitsWithoutPadding)/8, 8);
    txtBitsDecimal = bi2de(txtBits, 'left-msb');
    msg = char(txtBitsDecimal); % Convert ASCII to text
end

%% Plot BER vs SNR
figure(1);
semilogy(SNR, BER, '--');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('SNR vs BER for 8-PSK Modulation');
grid on;