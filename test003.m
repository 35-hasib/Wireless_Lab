%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amplitude Mod - DeMod

clc; close all; clear all;
fs = 1000; fm = 10; fc = 100;
m = 0.8;
t = 0 : 1/fs : 1;
m_s = sin(2*pi*fm*t);
c_s = cos(2*pi*fc*t);
mod_s = (1 + m*m_s).*c_s;

%% noise generation
%snr = 1;
%ps = mean(mod_s.^2);
%pn = ps/(10^(snr/10));
%noise = sqrt(pn/2) * randn(size(mod_s));
%n_mod_s = mod_s + noise;
n_mod_s = noisy_signal(mod_s);

%% Demodulation
demod_s = 2*n_mod_s .* c_s;
demod_s = lowpass(demod_s,fm,fs);

%% Figure
figure
plot(m_s)
figure
plot(c_s)
figure
plot(mod_s)
figure
plot(n_mod_s)
figure
plot(demod_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% BPSK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all;
bit = 5;
bit_s = abs(round(rand(1,bit)))
%%
fs = 100;
t = 0 : 1/fs: 1;
Eb = 1;
Tb = 1;
fc = 4;

%%
c_s = cos(2*pi*fc*t);
subplot(311)
plot(c_s)

%% PSK modulation
ttt = sqrt(2*Eb/Tb) * c_s;
TX = [];
for m = 1 : bit
    TX = [TX (-1)^bit_s(m) * (-1) *ttt];
end
%% Figures
subplot(312)
plot(TX)

RX = noisy_signal(TX);
subplot(313)
plot(RX)

%%
LO = sqrt(2/Tb) .* c_s;
det_bit_s = [];
cs = [];
for i = 1:bit
    temp = RX([(i-1)*101+1 : (i-1)*101+101]);
    s = sum(temp .*LO);
    cs = [cs s];
    if s>0
        det_bit_s = [det_bit_s 1];
    else
        det_bit_s = [det_bit_s 0];
    end
end

figure
subplot(311)
stem(cs)

subplot(312)
stem(bit_s)

subplot(313)
stem(det_bit_s)

bit_error = sum(abs(det_bit_s-bit_s))
bit_error_rate = bit_error/bit * 100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SNR vs BER %%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Importent Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Custom Functions
function result = de2bi(varargin)
    % Custom de2bi function
    if nargin == 2
        num = varargin{1};
        bits = varargin{2};
        result = zeros(length(num), bits);
        for i = 1:length(num)
            for j = 1:bits
                result(i, bits - j + 1) = mod(num(i), 2);
                num(i) = floor(num(i) / 2);
            end
        end
    else
        error('Custom de2bi function requires 2 inputs.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = bi2de(varargin)
    % Custom bi2de function
    if nargin == 2
        bits = varargin{1};
        result = zeros(size(bits, 1), 1);
        for i = 1:size(bits, 1)
            for j = 1:size(bits, 2)
                result(i) = result(i) + bits(i, j) * 2^(size(bits, 2) - j);
            end
        end
    else
        error('Custom bi2de function requires 2 inputs.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = pskmod(data, M)
    % Custom pskmod function
    phase = 2 * pi * data / M;
    result = cos(phase) + 1i * sin(phase);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = pskdemod(data, M)
    % Custom pskdemod function
    phase = angle(data);
    result = round(mod(phase * M / (2 * pi), M));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = awgn(data, snr, ~)
    % Custom awgn function
    signal_power = mean(abs(data).^2);
    noise_power = signal_power / (10^(snr / 10));
    noise = sqrt(noise_power / 2) * (randn(size(data)) + 1i * randn(size(data)));
    result = data + noise;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [noe, ber] = biterr(data1, data2)
    % Custom biterr function
    noe = sum(data1 ~= data2);
    ber = noe / length(data1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











