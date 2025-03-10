function f = noisy_signal(x)
    snr = 1;
    ps = mean(x .^2);
    pn = ps/(10^(snr/10));
    noise = sqrt(pn/2) * randn(size(x));
    f = x + noise;
end