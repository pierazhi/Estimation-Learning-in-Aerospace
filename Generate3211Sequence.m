clc; clearvars; close all;

%% ------------------------------------------------------------------------
% 3-2-1-1 Input Generator for Aircraft System Identification
% ------------------------------------------------------------------------
% This script generates a classical 3-2-1-1 excitation maneuver:
% - Pulse widths follow a 3:2:1:1 time ratio
% - The "2" pulse duration = half the period of expected dominant mode
% - Output: time vector and amplitude signal (repeated multiple times)
% ------------------------------------------------------------------------

% Sampling time
Ts = 0.0040;  % [s]

% Target frequency (dominant mode) in rad/s
omega_n = 3.2405;  % [rad/s]

% Duration of the "2" pulse = half-period of dominant mode
durstep = pi / omega_n;  % ~0.97s

% Define amplitude pattern for 3-2-1-1 sequence
amplitudes = [0,  0.1, -0.1,  0.1, -0.1,  0,  0];  % [neutral, 3, 2, 1, 1, neutral, neutral]

% Define time duration for each pulse
% Pulse width ratio: 3-2-1-1 → scaled using durstep
time_pattern = [1, 3*durstep, 2*durstep, durstep, durstep, 1, 1];  % [s]

% Total cycle duration
cycle_duration = sum(time_pattern);  % [s]

% Number of repetitions
num_repeats = 10;

% Convert one cycle to discrete samples
num_samples_per_cycle = round(cycle_duration / Ts);
time_single = (0:num_samples_per_cycle-1) * Ts;
signal_single = zeros(1, num_samples_per_cycle);

% Generate one cycle of the maneuver
start_idx = 1;
for i = 1:length(time_pattern)
    duration_samples = round(time_pattern(i) / Ts);
    end_idx = start_idx + duration_samples - 1;
    
    if end_idx > length(signal_single)
        end_idx = length(signal_single);
    end
    
    signal_single(start_idx:end_idx) = amplitudes(i);
    start_idx = end_idx + 1;
end

% Repeat full cycle
signal = repmat(signal_single, 1, num_repeats);
time = (0:length(signal)-1) * Ts;

% Plot 3-2-1-1 maneuver
figure;
stairs(time, signal, 'LineWidth', 1.5);
title('3-2-1-1 Maneuver (Normalized Pitching Moment)');
xlabel('Time [s]');
ylabel('Amplitude');
grid on;

% Save result
Excitation3211 = [time', signal'];
save('Excitation3211.mat', 'Excitation3211');

%% ------------------------------------------------------------------------
% Optional: Power Spectrum Analysis
% ------------------------------------------------------------------------

[Pxx, F] = periodogram(signal, [], [], 1/Ts);
figure;
plot(F, 10*log10(Pxx), 'LineWidth', 1.5);
title('Power Spectrum of 3-2-1-1 Input');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
xlim([0, 10]);
grid on;
