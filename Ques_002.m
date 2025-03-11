clc;
clear;
close all;

% Binary sequence
binary_sequence = [1 0 1 1 0 0 1 0 0 1];
number_bits = length(binary_sequence);

% Parameters
Eb = 1; % Energy per bit
Tb = 1; % Bit duration
Ac = 1; % Carrier amplitude
fc = 4 / Tb; % Carrier frequency
tf = 99; % Time factor
t = 0:1/tf:1; % Time vector for one bit
tn = 0:1/(tf+1):number_bits; % Time vector for all bits
tt = tn(1, 2:end); % Adjusted time vector

% Carrier signal
wc = 2 * pi * fc;
xc = Ac * cos(wc * t);

% Polar NRZ signal generation
NRZ = [];
for m = 1:number_bits
    if binary_sequence(m) == 1
        NRZ = [NRZ ones(1, length(t))];
    else
        NRZ = [NRZ -ones(1, length(t))];
    end
end

% BPSK signal generation
TX = [];
for n = 1:number_bits
    if binary_sequence(n) == 1
        TX = [TX sqrt(2 * Eb / Tb) * cos(2 * pi * fc * t)];
    else
        TX = [TX -sqrt(2 * Eb / Tb) * cos(2 * pi * fc * t)];
    end
end

% Add AWGN noise
SNR = 1; % Signal-to-Noise Ratio in dB
Ps = mean(abs(TX).^2); % Signal power
Pn = Ps / (10^(SNR / 10)); % Noise power
noise = sqrt(Pn) * randn(1, length(TX)); % Generate noise
RX = TX + noise; % Add noise to the signal

% Coherent demodulation
LO = sqrt(2 / Tb) * cos(2 * pi * fc * t); % Local oscillator
BINSEQDET = [];
CS = [];
for n = 1:number_bits
    temp = RX((n - 1) * (tf + 1) + 1:n * (tf + 1));
    S = sum(temp .* LO); % Correlation
    CS = [CS S];
    if S > 0
        BINSEQDET = [BINSEQDET 1];
    else
        BINSEQDET = [BINSEQDET 0];
    end
end

% Calculate bit errors
bit_errors = sum(abs(BINSEQDET - binary_sequence));
disp(['Total Bit Errors: ', num2str(bit_errors)]);

% Plotting
figure(1);
subplot(2, 1, 1);
plot(tt, NRZ);
title('Polar NRZ Input Binary Sequence');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(tt, TX(1, 1:length(tt)));
title('BPSK Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');

figure(2);
subplot(2, 1, 1);
plot(tt, RX(1, 1:length(tt)));
title('Received BPSK Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2, 1, 2);
stem(CS);
title('Output of the Coherent Correlation Receiver');

figure(3);
subplot(2, 1, 1);
stem(binary_sequence, 'Linewidth', 2);
title('Input Binary Sequence');
subplot(2, 1, 2);
stem(BINSEQDET, 'Linewidth', 2);
title('Detected Binary Sequence');