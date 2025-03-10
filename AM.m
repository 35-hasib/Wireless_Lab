
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