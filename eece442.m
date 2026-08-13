% EECE 442 Project 
% Runs the 4 modules in order. sample(), reconstruct() and quan() 
clear; clc; close all;

%% 1 Sampling
fprintf('1) Sampling\n');

f0 = 100;  fN = 2*f0;% test tone and its Nyquist rate
T  = 0.05;  fdense = 100*f0;% window length and the fine "analog" grid
t  = 0:1/fdense:T;
xt = cos(2*pi*f0*t);

fs_low  = 0.5*fN; % under Nyquist, should alias
fs_high = 2*fN; % over Nyquist, should be fine

[ts_lo, xs_lo] = sample(t, xt, fs_low);
[ts_hi, xs_hi] = sample(t, xt, fs_high);
xr_lo = reconstruct(t, xs_lo, fs_low);
xr_hi = reconstruct(t, xs_hi, fs_high);

figure('Name','M1: cosine sampling');
subplot(2,1,1);
plot(t, xt, 'k'); hold on;
stem(ts_lo, xs_lo, 'r', 'filled');
plot(t, xr_lo, 'r--', 'LineWidth', 1.2);
title(sprintf('fs = 0.5 f_N = %g Hz (below Nyquist -> aliasing)', fs_low));
legend('original','samples','reconstruction'); xlabel('t [s]'); grid on;
subplot(2,1,2);
plot(t, xt, 'k'); hold on;
stem(ts_hi, xs_hi, 'b', 'filled');
plot(t, xr_hi, 'b--', 'LineWidth', 1.2);
title(sprintf('fs = 2 f_N = %g Hz (above Nyquist -> proper sampling)', fs_high));
legend('original','samples','reconstruction'); xlabel('t [s]'); grid on;

in = t > 0.1*T & t < 0.9*T;% skip the edges, the sinc sum is off there
fprintf('Cosine MSE (fs = 0.5 fN): %.4e\n', mean((xt(in)-xr_lo(in)).^2));
fprintf('Cosine MSE (fs = 2.0 fN): %.4e\n', mean((xt(in)-xr_hi(in)).^2));

% same aliasing story but in frequency: the low-rate one should peak at 0 Hz
Nf = 2^14;  f = (0:Nf-1)*fdense/Nf;
Xhi = abs(fft(xr_hi, Nf));  Xlo = abs(fft(xr_lo, Nf));
figure('Name','M1: aliasing in frequency');
plot(f, Xhi/max(Xhi), 'b', 'LineWidth', 1.1); hold on;
plot(f, Xlo/max(Xlo), 'r--', 'LineWidth', 1.1);
xlim([0 300]); grid on; xlabel('f [Hz]'); ylabel('|X(f)| (normalized)');
legend('reconstructed at 2f_N: peak at 100 Hz', ...
       'reconstructed at 0.5f_N: alias at 0 Hz');
title('Aliasing in the frequency domain');

% message: two tones, like a phone keypad tone. Highest freq is 400 Hz.
f1 = 250;  f2 = 400;  fmax = max(f1,f2);
mN   = 2*fmax; % its Nyquist rate
m    = 0.7*cos(2*pi*f1*t) + 0.5*sin(2*pi*f2*t);
fs_m = 2*mN;% sample well above it

[tm, ms] = sample(t, m, fs_m);
mr = reconstruct(t, ms, fs_m);

figure('Name','M1: message signal');
plot(t, m, 'k'); hold on;
plot(t, mr, 'g--', 'LineWidth', 1.2);
stem(tm, ms, 'g', 'filled', 'MarkerSize', 3);
title(sprintf('Message m(t), fs = %g Hz (Nyquist rate = %g Hz)', fs_m, mN));
legend('m(t)','reconstruction','samples'); xlabel('t [s]'); grid on;
fprintf('Message reconstruction MSE: %.4e\n', mean((m(in)-mr(in)).^2));

% longer 2 s record so the stats in modules 2-3 settle down
t_long = 0:1/fdense:2;
m_long = 0.7*cos(2*pi*f1*t_long) + 0.5*sin(2*pi*f2*t_long);
[tm2, x] = sample(t_long, m_long, fs_m);

%% 2) Quantization
fprintf('\n2) Quantization\n');
xmin = min(x);  xmax = max(x);

% two levels: one threshold in the middle, levels at the two half-midpoints
thr2 = (xmin + xmax)/2;
lev2 = [ (xmin + thr2)/2 , (thr2 + xmax)/2 ];
xq2  = quan(x, thr2, lev2);

figure('Name','M2: two-level');
plot(tm2, x, 'k.-'); hold on;
stairs(tm2, xq2, 'r', 'LineWidth', 1.1);
title('Two-level quantization');
legend('sampled x','quantized x_q'); xlabel('t [s]'); xlim([0 0.05]); grid on;

Ms  = [4 16 32]; % try low and high rates
mse = zeros(size(Ms));
fprintf('%-6s %-12s\n','M','MSE');
for i = 1:length(Ms)
    M = Ms(i);
    D = (xmax - xmin)/M; % step size
    thr = xmin + D*(1:M-1);% thresholds sit on the interval edges
    lev = xmin + D*((1:M) - 0.5);% levels sit at the interval midpoints
    xq  = quan(x, thr, lev);
    mse(i) = mean((x - xq).^2);

    fprintf('%-6d %-12.4e\n', M, mse(i));
    fprintf('  thresholds: %s\n', mat2str(round(thr,4)));
    fprintf('  levels    : %s\n', mat2str(round(lev,4)));

    figure('Name',sprintf('M2: M = %d', M));
    plot(tm2, x, 'k.-'); hold on;
    stairs(tm2, xq, 'b', 'LineWidth', 1.1);
    title(sprintf('Uniform quantization, M = %d, MSE = %.3e', M, mse(i)));
    legend('sampled x','quantized x_q'); xlabel('t [s]'); xlim([0 0.05]); grid on;
end

figure('Name','M2: MSE vs M');
semilogy(Ms, mse, 'o-'); grid on;
xlabel('Number of levels M'); ylabel('MSE');
title('MSE vs quantization rate');

sqnr = 10*log10(mean(x.^2)./mse); % should go up ~6 dB every time we add a bit
fprintf('SQNR [dB] for M = [4 16 32]: %s\n', mat2str(round(sqnr,2)));

% build the M=16 symbol stream (which bin each sample lands in, 1..16)
M16 = 16;  D = (xmax - xmin)/M16;
thr16    = xmin + D*(1:M16-1);
lev_init = xmin + D*((1:M16) - 0.5);
stream = ones(size(x));
for k = 1:length(thr16), stream = stream + (x > thr16(k)); end
stream = stream(:);

% Lloyd-Max: same M=16 but levels moved to minimize MSE for this signal.
% Start from the uniform levels and let lloyds() refine them, then quantize
% with the same quan() function to show it handles any threshold set.
[thr_lm, lev_lm] = lloyds(x, lev_init);
xq_lm  = quan(x, thr_lm, lev_lm);
mse_lm = mean((x - xq_lm).^2);
fprintf('Lloyd-Max M = 16:  MSE = %.4e  (uniform M = 16: %.4e)\n', mse_lm, mse(2));
fprintf('  LM thresholds: %s\n', mat2str(round(thr_lm,4)));
fprintf('  LM levels    : %s\n', mat2str(round(lev_lm,4)));
sym_lm = ones(size(x));
for k = 1:length(thr_lm), sym_lm = sym_lm + (x > thr_lm(k)); end
sym_lm = sym_lm(:);

% overlay the signal with both M=16 versions so the difference is visible
figure('Name','M2: chain overlay');
plot(tm2, x, 'k', 'LineWidth', 1.0); hold on;
stairs(tm2, quan(x, thr16, lev_init), 'b', 'LineWidth', 1.0);
stairs(tm2, xq_lm, 'r', 'LineWidth', 1.0);
xlim([0 0.01]); grid on; xlabel('t [s]');
legend('sampled m(t)','uniform M = 16','Lloyd-Max M = 16');
title('Message vs uniform and Lloyd-Max quantization (M = 16)');

%% 3) Source Coding
fprintf('\n3) Source Coding\n');
N = length(stream);

Lfix = ceil(log2(M16)); % fixed-length: same bits for every symbol
encFix = reshape(dec2bin(stream-1, 4).' - '0', [], 1); % actually encode it, 4 bits each
fprintf('(a) Fixed-length baseline: %d bits/symbol\n', Lfix);

A = unique(stream).';
p = zeros(size(A));
for i = 1:length(A)
    p(i) = sum(stream == A(i)) / N;% how often each symbol shows up
end
H = -sum(p .* log2(p)); % entropy = best possible bits/symbol
fprintf('(b) Entropy H(A) = %.4f bits/symbol\n', H);

figure('Name','M3: empirical distribution');
bar(A, p); xlabel('symbol index a'); ylabel('p(a)');
title(sprintf('Empirical distribution, H(A) = %.3f bits', H)); grid on;

dict = huffmandict(A, p); % build the Huffman code from those frequencies
enc  = huffmanenco(stream, dict);
fprintf('    Fixed-length stream: %d bits; Huffman stream: %d bits\n', ...
        numel(encFix), length(enc));
dec  = huffmandeco(enc, dict);
lossless = isequal(dec(:), stream(:)); % decoded should match exactly
fprintf('(c) Lossless reconstruction: %s\n', string(lossless));
assert(lossless, 'Huffman decode does not match source stream.');

codelens = cellfun(@length, dict(:,2)).';
Lbar     = sum(p .* codelens); % average bits/symbol Huffman actually used
fprintf('(d) Lbar = %.4f bits/symbol\n', Lbar);
fprintf('    H(A) <= Lbar < H(A)+1: %.4f <= %.4f < %.4f\n', H, Lbar, H+1);
fprintf('    Compression vs fixed-length: %.1f%%\n', 100*Lbar/Lfix);

% same coding on the Lloyd-Max stream, to compare the two quantizers end to end
A2 = unique(sym_lm).';
p2 = zeros(size(A2));
for i = 1:length(A2)
    p2(i) = sum(sym_lm == A2(i)) / N;
end
H2 = -sum(p2 .* log2(p2));
dict2 = huffmandict(A2, p2);
enc2  = huffmanenco(sym_lm, dict2);
dec2  = huffmandeco(enc2, dict2);
assert(isequal(dec2(:), sym_lm(:)), 'LM Huffman decode mismatch.');
Lbar2 = sum(p2 .* cellfun(@length, dict2(:,2)).');
fprintf('Lloyd-Max stream:  H = %.4f, Lbar = %.4f bits/symbol, Huffman = %d bits\n', ...
        H2, Lbar2, length(enc2));

%% 4) BPSK/QPSK over AWGN
fprintf('\n4) BPSK/QPSK over AWGN\n');
rng(442); % fixed seed so results repeat
Nb   = 1e6;
bits = randi([0 1], Nb, 1);

EbN0_dB = 0:1:10;
berBPSK = zeros(size(EbN0_dB));
berQPSK = zeros(size(EbN0_dB));

for i = 1:length(EbN0_dB)
    % BPSK: 1 bit per symbol, so symbol SNR = Eb/N0
    xB = pskmod(bits, 2);
    yB = awgn(xB, EbN0_dB(i), 'measured');
    bB = pskdemod(yB, 2);
    berBPSK(i) = sum(bB ~= bits) / Nb;

    % QPSK: 2 bits per symbol, so bump the SNR by 10log10(2) to match per-bit energy
    symsTx = bit2int(bits, 2); % pair up bits into 0..3
    xQ = pskmod(symsTx, 4, pi/4, 'gray');
    yQ = awgn(xQ, EbN0_dB(i) + 10*log10(2), 'measured');
    symsRx = pskdemod(yQ, 4, pi/4, 'gray');
    bQ = int2bit(symsRx, 2);% back to bits
    berQPSK(i) = sum(bQ ~= bits) / Nb;
end

EbN0  = 10.^(EbN0_dB/10);
berTh = qfunc(sqrt(2*EbN0));% theory curve, same for both here

figure('Name','P2: BER vs Eb/N0');
semilogy(EbN0_dB, berBPSK, 'bo-'); hold on;
semilogy(EbN0_dB, berQPSK, 'rs-');
semilogy(EbN0_dB, berTh,   'k--', 'LineWidth', 1.2);
grid on; xlabel('E_b/N_0 [dB]'); ylabel('BER');
legend('BPSK (sim)','QPSK (sim)','Q(\surd(2E_b/N_0)) theory');
title('BER estimation over AWGN, 10^6 bits');

fprintf('%-10s %-12s %-12s %-12s\n','Eb/N0 dB','BPSK sim','QPSK sim','Theory');
fprintf('%-10d %-12.3e %-12.3e %-12.3e\n', [EbN0_dB; berBPSK; berQPSK; berTh]);

% show the received QPSK points at low vs high SNR: the clouds tighten up
figure('Name','P2: QPSK constellations');
snrs = [2 10];
for s = 1:2
    yQs = awgn(pskmod(bit2int(bits(1:2e4), 2), 4, pi/4, 'gray'), ...
               snrs(s) + 10*log10(2), 'measured');
    subplot(1,2,s);
    plot(real(yQs), imag(yQs), '.', 'MarkerSize', 3); axis square; grid on;
    axis([-2.5 2.5 -2.5 2.5]);% same axes on both so the shrink is honest
    xlabel('I'); ylabel('Q');
    title(sprintf('QPSK received, E_b/N_0 = %d dB', snrs(s)));
end

fprintf('\nAll modules completed.\n');

%% functions
function [t_sample, x_sample] = sample(t, xt, fs)
% grab samples at 1/fs spacing, reading values off the fine grid
Ts = 1/fs;
t_sample = t(1):Ts:t(end);
x_sample = interp1(t, xt, t_sample, 'linear');
end

function xr = reconstruct(t, x_sample, fs)
% rebuild the signal with the sinc formula: sum of shifted sincs
Ts = 1/fs;
n  = 0:length(x_sample)-1;
tn = t(1) + n*Ts;% where each sample sits in time
xr = x_sample * sinc(fs * (t - tn.'));
end

function xq = quan(x, thr, levels)
% send each sample to a level: count how many thresholds it passed
thr = sort(thr(:)).';
idx = ones(size(x));
for k = 1:length(thr)
    idx = idx + (x > thr(k));
end
xq = levels(idx);
end