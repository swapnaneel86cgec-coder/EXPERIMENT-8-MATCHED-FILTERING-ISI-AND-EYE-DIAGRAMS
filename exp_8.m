% SAVE THIS FILE EXACTLY AS: experiment_8.m
% Do not use spaces, hyphens, or start the filename with a number.

% Experiment 8 - Matched filtering, ISI and eye diagrams
clear; clc; close all;

% ==========================================
% 1. System Parameters
% ==========================================
num_symbols = 5000;
sps = 8;          % Samples per symbol
span = 6;         % Filter span in symbols
beta = 0.3;       % Roll-off factor
EbNo_dB = 10;     % SNR for the basic time-domain plots

% ==========================================
% 2. Baseband Tx & Rx 
% ==========================================
% Generate BPSK symbols (+1, -1)
bits = randi([0, 1], 1, num_symbols);
symbols = 2 * bits - 1;

% Upsample
tx_upsampled = zeros(1, num_symbols * sps);
tx_upsampled(1:sps:end) = symbols;

% Pulse Shaping (Tx)
h_rrc = custom_rrc(sps, span, beta);
tx_signal = conv(tx_upsampled, h_rrc, 'full');

% AWGN Channel
signal_power = mean(tx_signal.^2);
snr_linear = 10^(EbNo_dB / 10);
noise_variance = signal_power / (2 * snr_linear);
noise = sqrt(noise_variance) * randn(1, length(tx_signal));
rx_signal = tx_signal + noise;

% Matched Filter (Rx)
rx_matched = conv(rx_signal, h_rrc, 'full');

% Delay Compensation
filter_delay = span * sps; % Total delay from Tx and Rx filters
rx_compensated = rx_matched(filter_delay + 1 : filter_delay + (num_symbols * sps));

disp('--- VALIDATION ---');
fprintf('Total Cascade Delay computed as: %d samples\n', filter_delay);

% ==========================================
% 3. Required Visualizations 
% ==========================================

% Viz 1: Transmit / Matched-filter output comparison
figure('Name', 'Tx vs Rx Output', 'Position', [100, 100, 800, 400]);
plot(1:100, tx_upsampled(1:100), 'o'); hold on;
plot(1:100, rx_compensated(1:100), 'LineWidth', 1.5);
title('Transmit Impulses vs. Matched Filter Output');
xlabel('Samples'); ylabel('Amplitude');
legend('Upsampled Impulses (Tx)', 'Matched Filter Output');
grid on; hold off;

% Viz 2: Eye Diagrams (Ideal vs Multipath)
plot_eye_diagram(rx_compensated, sps, sprintf('Eye Diagram (AWGN, SNR=%ddB)', EbNo_dB));

% Introduce Multipath ISI
multipath_channel = [1.0, 0.4]; 
rx_multipath = conv(tx_signal, multipath_channel, 'full');
rx_multipath = rx_multipath(1:length(tx_signal)) + noise;
rx_matched_isi = conv(rx_multipath, h_rrc, 'full');
rx_isi_comp = rx_matched_isi(filter_delay + 1 : filter_delay + (num_symbols * sps));
plot_eye_diagram(rx_isi_comp, sps, 'Eye Diagram with Multipath ISI');

% Viz 3: Eye height versus SNR
snr_range = 0:2:14;
eye_heights = zeros(1, length(snr_range));

for k = 1:length(snr_range)
    snr = snr_range(k);
    n_var = signal_power / (2 * (10^(snr / 10)));
    n = sqrt(n_var) * randn(1, length(tx_signal));
    r = tx_signal + n;
    rm = conv(r, h_rrc, 'full');
    rm_comp = rm(filter_delay + 1 : filter_delay + (num_symbols * sps));
    
    opt_samples = rm_comp(1:sps:end);
    eye_heights(k) = 2 * mean(abs(opt_samples)) - 6 * std(opt_samples); 
end

figure('Name', 'Eye Height vs SNR', 'Position', [100, 100, 700, 400]);
plot(snr_range, eye_heights, '-o', 'LineWidth', 1.5);
title('Estimated Eye Height vs SNR');
xlabel('SNR (dB)'); ylabel('Eye Height');
grid on;

% Viz 4: BER versus timing offset
timing_offsets = 0:(sps-1);
ber_array = zeros(1, length(timing_offsets));

for k = 1:length(timing_offsets)
    offset = timing_offsets(k);
    sampled_rx = rx_compensated((offset + 1):sps:end);
    sampled_rx = sampled_rx(1:num_symbols);
    detected_bits = double(sampled_rx > 0);
    ber_array(k) = sum(bits ~= detected_bits) / num_symbols;
end

figure('Name', 'BER vs Timing Offset', 'Position', [100, 100, 700, 400]);
semilogy(timing_offsets, ber_array, '-sr', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
title('BER vs Timing Offset (Sampling Instant)');
xlabel('Timing Offset (Samples from Optimum)'); ylabel('Bit Error Rate (BER)');
grid on;


% ==========================================
% 4. Helper Functions
% ==========================================

function h = custom_rrc(sps, span, beta)
    t = -span/2 : 1/sps : span/2;
    t = t + 1e-8; 
    h = zeros(1, length(t));
    
    for i = 1:length(t)
        tc = t(i);
        if abs(tc) < 1e-7
            h(i) = 1.0 - beta + (4 * beta / pi);
        elseif beta ~= 0 && abs(abs(tc) - 1/(4*beta)) < 1e-7
            h(i) = (beta / sqrt(2)) * (((1 + 2/pi) * sin(pi/(4*beta))) + ((1 - 2/pi) * cos(pi/(4*beta))));
        else
            num = sin(pi*tc*(1-beta)) + 4*beta*tc.*cos(pi*tc*(1+beta));
            den = pi*tc * (1 - (4*beta*tc)^2);
            h(i) = num / den;
        end
    end
    h = h / sqrt(sum(h.^2)); 
end

function plot_eye_diagram(signal, sps, fig_title)
    samples_per_trace = 2 * sps; 
    num_traces = floor(length(signal) / samples_per_trace);
    
    figure('Name', fig_title, 'Position', [100, 100, 700, 400]);
    hold on;
    time_axis = linspace(-1, 1, samples_per_trace);
    traces_to_plot = min(num_traces, 500);
    
    for i = 0:(traces_to_plot - 1)
        start_idx = i * sps + 1;
        end_idx = start_idx + samples_per_trace - 1;
        
        if end_idx <= length(signal)
            plot(time_axis, signal(start_idx:end_idx), 'Color', [0 0 1 0.1]);
        end
    end
    
    title(fig_title);
    xlabel('Symbol Intervals (T)'); ylabel('Amplitude');
    grid on; hold off;
end