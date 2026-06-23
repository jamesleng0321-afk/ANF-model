%% ========================================================================
% Evaluate_Facilitation_JoshiFig6a.m
% Reproduction of Joshi 2017 Fig 6a: Responses to subthreshold paired pulse
% Evaluates probe threshold shift as a function of inter-pulse delay.
% Comparison: Joshi 2017 vs Bruce 2024
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Global Parameters
% ------------------------------------------------------------------------
Fs = 1e6;
NoiseAlpha = 0.8;
nTrials = 25; % Monte Carlo trials for threshold estimation

% 100 us/phase biphasic pulse (cathodic first, matches Dynes 1996 model approximations)
phase_len_us = 100;
SinglePulse = [0, -1*ones(1, phase_len_us), 1*ones(1, phase_len_us), 0];
pulse_dur_s = length(SinglePulse) / Fs;

% Conditioner levels according to Joshi 2017 Fig 6a: -0.9 dB and -2.0 dB
cond_levels_dB = [-0.9, -2.0];

% Inter-pulse delays (gap between offset of conditioner and onset of probe)
% Sweep from 50 us to 4000 us
delays_us = [50, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000];
delays_s  = delays_us * 1e-6;

%% ------------------------------------------------------------------------
% 2. Calculate Baseline Single-Pulse Threshold (I50)
% ------------------------------------------------------------------------
fprintf('=== Phase 1: Calculating Baseline I50 Thresholds ===\n');

% Range for binary search [0.2 mA to 2.0 mA]
min_I = 0.2e-3; max_I = 2.0e-3;

I50_Joshi = fast_threshold_search(@Model_PulseTrain, 0, 0, SinglePulse, Fs, NoiseAlpha, nTrials, min_I, max_I);
I50_Bruce = fast_threshold_search(@Model_PulseTrain_Bruce, 0, 0, SinglePulse, Fs, NoiseAlpha, nTrials, min_I, max_I);

fprintf('  -> Joshi 2017 Baseline I50: %.3f mA\n', I50_Joshi * 1000);
fprintf('  -> Bruce 2024 Baseline I50: %.3f mA\n', I50_Bruce * 1000);

%% ------------------------------------------------------------------------
% 3. Main Simulation Loop: Paired Pulse Facilitation
% ------------------------------------------------------------------------
fprintf('\n=== Phase 2: Simulating Paired Pulse Thresholds ===\n');

% Preallocate results: 
% Dimensions: (Model: 1=Joshi, 2=Bruce) x (Conditioner Level) x (Delays)
threshold_shifts_dB = zeros(2, length(cond_levels_dB), length(delays_us));

for lvl_idx = 1:length(cond_levels_dB)
    cond_dB = cond_levels_dB(lvl_idx);
    fprintf('\nProcessing Conditioner Level: %.1f dB...\n', cond_dB);
    
    % Calculate exact amplitude for the conditioner pulse
    amp_cond_J = I50_Joshi * 10^(cond_dB / 20);
    amp_cond_B = I50_Bruce * 10^(cond_dB / 20);
    
    for d_idx = 1:length(delays_us)
        gap_s = delays_s(d_idx);
        
        % 1. Joshi 2017 Probe Threshold
        probe_I50_J = fast_threshold_search(@Model_PulseTrain, amp_cond_J, gap_s, SinglePulse, Fs, NoiseAlpha, nTrials, 0.1e-3, 1.5e-3);
        threshold_shifts_dB(1, lvl_idx, d_idx) = 20 * log10(probe_I50_J / I50_Joshi);
        
        % 2. Bruce 2024 Probe Threshold
        probe_I50_B = fast_threshold_search(@Model_PulseTrain_Bruce, amp_cond_B, gap_s, SinglePulse, Fs, NoiseAlpha, nTrials, 0.1e-3, 1.5e-3);
        threshold_shifts_dB(2, lvl_idx, d_idx) = 20 * log10(probe_I50_B / I50_Bruce);
    end
end

%% ------------------------------------------------------------------------
% 4. Plotting (Recreating Joshi 2017 Fig 6a Style)
% ------------------------------------------------------------------------
fprintf('\n=== Phase 3: Plotting Results ===\n');

fig = figure("Name", "Subthreshold Paired Pulse Stimulation", "DefaultAxesFontSize", 13);
fig.Position(3:4) = [800 600];
hold on;

% Plot reference line at 0 dB (No threshold shift)
plot([0 4000], [0 0], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Define visual styles
colors = ["#D95319", "#0072BD"]; % Orange for -0.9dB, Blue for -2.0dB
markers_Joshi = ['o', 'o'];
markers_Bruce = ['s', 's'];

lgh = [];
lg_labels = {};

for lvl_idx = 1:length(cond_levels_dB)
    % Plot Joshi 2017 (Solid lines with Circles)
    l_j = plot(delays_us, squeeze(threshold_shifts_dB(1, lvl_idx, :)), ...
        '-', 'Marker', markers_Joshi(lvl_idx), 'LineWidth', 2, 'MarkerSize', 7, ...
        'Color', colors(lvl_idx), 'MarkerFaceColor', colors(lvl_idx));
    
    % Plot Bruce 2024 (Dashed lines with Squares)
    l_b = plot(delays_us, squeeze(threshold_shifts_dB(2, lvl_idx, :)), ...
        '--', 'Marker', markers_Bruce(lvl_idx), 'LineWidth', 2, 'MarkerSize', 7, ...
        'Color', colors(lvl_idx), 'MarkerFaceColor', 'none');
        
    lgh = [lgh, l_j, l_b];
    lg_labels = [lg_labels, ...
        sprintf('%.1f dB (Joshi 2017)', cond_levels_dB(lvl_idx)), ...
        sprintf('%.1f dB (Bruce 2024)', cond_levels_dB(lvl_idx))];
end

% Formatting to match the original paper
grid on; box off;
xlim([0 4000]);
ylim([-14 2]);

xlabel('Inter-pulse delay (\mus)', 'FontWeight', 'bold');
ylabel('Threshold shift (dB)', 'FontWeight', 'bold');
title('Probe Threshold Shift vs Inter-pulse Delay');

legend(lgh, lg_labels, 'Location', 'southeast', 'FontSize', 11);

disp('Simulation complete.');

%% ========================================================================
% Helper Functions
% ========================================================================

function threshold = fast_threshold_search(model_func, cond_amp, gap_s, pulse, Fs, NoiseAlpha, nTrials, min_I, max_I)
    % Performs a binary search to find the amplitude that yields 50% spiking probability
    target_prob = 0.5;
    tolerance = 0.05;
    max_iters = 8; % 8 iterations is usually enough for 0.01 mA precision
    
    pulse_dur_s = length(pulse) / Fs;
    
    for iter = 1:max_iters
        test_amp = (min_I + max_I) / 2;
        
        % Build stimulus vector (Conditioner + Gap + Probe)
        % If cond_amp == 0, it behaves like finding single pulse threshold
        if cond_amp > 0
            gap_samples = round(gap_s * Fs);
            stim = [pulse * cond_amp, zeros(1, gap_samples), pulse * test_amp];
            probe_start_s = pulse_dur_s + gap_s;
        else
            stim = pulse * test_amp;
            probe_start_s = 0;
        end
        
        % Add padding at the end to allow spikes to occur
        stim = [stim, zeros(1, round(5e-3 * Fs))]; 
        
        spikes_cnt = 0;
        for tr = 1:nTrials
            p_n = Library.oneonfnoise(length(stim), NoiseAlpha);
            c_n = Library.oneonfnoise(length(stim), NoiseAlpha);
            [~, SpT, ~, ~, ~, ~] = model_func(stim, p_n, c_n, Fs);
            
            % Check if a spike occurred AFTER the probe onset
            if any(SpT >= probe_start_s)
                spikes_cnt = spikes_cnt + 1;
            end
        end
        
        current_prob = spikes_cnt / nTrials;
        
        % Binary search update rule
        if current_prob > target_prob
            max_I = test_amp; % Need less current
        else
            min_I = test_amp; % Need more current
        end
    end
    
    threshold = (min_I + max_I) / 2;
end