
clc; close all; clear all;
bit = 30;
bit_s = abs(round(rand(1,bit)));
%%
t = 0 : 1/100: 1;
Eb = 1;
Tb = 1;
fc = 4;

%%
c_s = cos(2*pi*fc*t);
figure
plot(c_s)
%% PSK modulation
ttt = sqrt(2*Eb/Tb) * c_s;
TX = [];
for m = 1 : bit
    if bit_s(m) == 1
        TX = [TX ttt];
    else
        TX = [TX -1*ttt];
    end
end

figure
plot(TX)

n_TX = noisy_signal(TX);
figure
plot(n_TX)


