%% ========================================================================
% Miller2011_PSTH_FE_Calibrated.m
% Miller 2011 Fig 1 / Boulet 2017 Fig 10 - PSTH Masking Recovery
% Evaluates Probe recovery after high-rate masker at specific Firing Efficiencies (FE)
% Comparison: Joshi 2017 vs Bruce 2024
% ========================================================================
clear; clc; close all;

%% 1. Parameters Configuration
Fs = 1e6;                       % Sampling frequency (1 MHz)
NoiseAlpha = 0.8;               % 1/f noise spectral shaping parameter
nTrials_Calib = 15;             % Monte Carlo trials for calibration phase
nTrials_Sim = 50;               % Monte Carlo trials for full simulation (increase for smoother PSTH)

% Masker settings: 200 ms duration, 5000 pps
masker_duration = 0.2;          % seconds
masker_rate     = 5000;         % pps
onset_window    = 12e-3;        % 12 ms onset window for FE calibration

% Probe settings: 250 ms duration, 100 pps
probe_duration  = 0.25;         % seconds
probe_rate      = 100;          % pps
probe_dB        = 0.0;          % Probe level at I50 (50% firing probability)

% Gap between Masker offset and Probe onset
gap_duration    = 0.2e-3;       % 0.2 ms

% Target Firing Efficiencies (FE) for the Masker
target_masker_FEs = [0.01, 0.50, 0.99, 0.9999]; 

% Pulse Definition: 40 us/phase symmetric biphasic, leading cathodic, 0 IPG
phase_len_us = 40;
pulse_ANF = [0, -1*ones(1, phase_len_us), 1*ones(1, phase_len_us), 0];

%% 2. Calculate Baseline Threshold (I50) for the Probe
disp('Calculating Single Pulse Baseline Threshold (I50) for Probe...');
% Assuming Library.FindThreshold and Library.FitNeuronDynamicRange are available
try
    [Level_curve, Prob_curve] = Library.FindThreshold([pulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
    [muSingle, ~] = Library.FitNeuronDynamicRange(Level_curve', Prob_curve);
    I50 = muSingle;
catch
    warning('Library functions not found. Using default I50 = 600 uA.');
    I50 = 600e-6;
end
fprintf('I50 found: %.4f uA\n', I50 * 1e6);

%% 3. Construct Base Stimuli Vectors
% Safe pulse generator to prevent array dimension mismatch
masker_base = build_safe_pulse_train(pulse_ANF, masker_rate, masker_duration, Fs);
probe_base  = build_safe_pulse_train(pulse_ANF, probe_rate, probe_duration, Fs);
gap_zeros   = zeros(1, round(gap_duration * Fs));

% Define timing boundaries for analysis
masker_end  = length(masker_base) / Fs;
probe_start = masker_end + gap_duration;
t_max       = probe_start + probe_duration;

%% 4. Auto-Calibration Phase (Binary Search for Masker FEs)
fprintf('\n--- Phase 1: Auto-Calibrating Masker Currents for Target FEs ---\n');

calibrated_amps_Joshi = zeros(1, length(target_masker_FEs));
calibrated_amps_Bruce = zeros(1, length(target_masker_FEs));

min_I = 0.5e-3; 
max_I = 3.0e-3; 

for fe_idx = 1:length(target_masker_FEs)
    target_FE = target_masker_FEs(fe_idx);
    target_rate = target_FE * masker_rate; 
    fprintf('\nTarget FE: %.2f%% (Onset Rate: %.1f spikes/s)\n', target_FE*100, target_rate);
    
    % Calibrate Joshi 2017
    calibrated_amps_Joshi(fe_idx) = calibrate_model_current(@Model_PulseTrain, ...
        pulse_ANF, masker_rate, onset_window, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials_Calib);
    fprintf('  -> Joshi 2017 calibrated current: %.3f mA\n', calibrated_amps_Joshi(fe_idx)*1000);
    
    % Calibrate Bruce 2024
    calibrated_amps_Bruce(fe_idx) = calibrate_model_current(@Model_PulseTrain_Bruce, ...
        pulse_ANF, masker_rate, onset_window, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials_Calib);
    fprintf('  -> Bruce 2024 calibrated current: %.3f mA\n', calibrated_amps_Bruce(fe_idx)*1000);
end

%% 5. Unmasked Baselines (Probe Only)
fprintf('\n--- Phase 2: Simulating Unmasked Probe Baselines ---\n');
amp_probe = I50 * 10^(probe_dB/20);
stim_unmasked = [zeros(1, length(masker_base) + length(gap_zeros)), probe_base * amp_probe];

[spikes_orig_unmasked, ~, sr_probe_orig_unm] = run_model_MC(@Model_PulseTrain, stim_unmasked, Fs, NoiseAlpha, nTrials_Sim, masker_end, probe_start, masker_duration, probe_duration);
[spikes_bruce_unmasked, ~, sr_probe_bruce_unm] = run_model_MC(@Model_PulseTrain_Bruce, stim_unmasked, Fs, NoiseAlpha, nTrials_Sim, masker_end, probe_start, masker_duration, probe_duration);

%% 6. Masked Simulation Loop
fprintf('\n--- Phase 3: Running Masked Accomodation Simulations ---\n');
n_levels = length(target_masker_FEs);

% Data structures to store results
res_orig  = struct('all_spikes', {}, 'mask_sr', {}, 'probe_sr', {});
res_bruce = struct('all_spikes', {}, 'mask_sr', {}, 'probe_sr', {});

for lvl_idx = 1:n_levels
    target_FE = target_masker_FEs(lvl_idx);
    fprintf('Processing Masker FE: %.2f%%\n', target_FE * 100);
    
    % Build masked stimuli using calibrated amplitudes
    stim_masked_joshi = [masker_base * calibrated_amps_Joshi(lvl_idx), gap_zeros, probe_base * amp_probe];
    stim_masked_bruce = [masker_base * calibrated_amps_Bruce(lvl_idx), gap_zeros, probe_base * amp_probe];
    
    % Run Original Model
    [sp, m_sr, p_sr] = run_model_MC(@Model_PulseTrain, stim_masked_joshi, Fs, NoiseAlpha, nTrials_Sim, masker_end, probe_start, masker_duration, probe_duration);
    res_orig(lvl_idx).all_spikes = sp; res_orig(lvl_idx).mask_sr = m_sr; res_orig(lvl_idx).probe_sr = p_sr;
    
    % Run Bruce Model
    [sp, m_sr, p_sr] = run_model_MC(@Model_PulseTrain_Bruce, stim_masked_bruce, Fs, NoiseAlpha, nTrials_Sim, masker_end, probe_start, masker_duration, probe_duration);
    res_bruce(lvl_idx).all_spikes = sp; res_bruce(lvl_idx).mask_sr = m_sr; res_bruce(lvl_idx).probe_sr = p_sr;
end

%% 7. Plotting PSTHs (Boulet 2017 Fig 10 Style)
% Configuration for PSTH binning
bin_width = 2e-3; % 2 ms bins for PSTH calculation
edges = 0 : bin_width : t_max;
t_centers = edges(1:end-1) + bin_width/2;

fig = figure("Name", "Miller2011_PSTH_FE_Comparison", "DefaultAxesFontSize", 11, "Color", "w");
fig.Position(3:4) = [1000, 900];
% 2 Columns: Original (Left) vs Bruce (Right), Rows: Unmasked + FEs
tiledlayout(n_levels + 1, 2, "TileSpacing", "compact", "Padding", "compact");

% Helper function for PSTH binning
get_psth = @(spikes) histcounts(spikes, edges) / (nTrials_Sim * bin_width);

% --- ROW 1: Unmasked Baselines ---
psth_orig_unm = get_psth(spikes_orig_unmasked);
psth_bruce_unm = get_psth(spikes_bruce_unmasked);
max_y = max([max(psth_orig_unm), max(psth_bruce_unm)]) * 1.2;

% Col 1: Original Unmasked
nexttile; hold on;
bar(t_centers, psth_orig_unm, 1, 'FaceColor', '#0072BD', 'EdgeColor', 'none');
ylabel(sprintf('Unmasked\nProbe SR: %.1f', sr_probe_orig_unm), 'FontWeight', 'bold');
title('Joshi 2017');
xlim([0, t_max]); ylim([0, max_y]); set(gca, 'XTickLabel', []); box off;

% Col 2: Bruce Unmasked
nexttile; hold on;
bar(t_centers, psth_bruce_unm, 1, 'FaceColor', '#D95319', 'EdgeColor', 'none');
ylabel(sprintf('Unmasked\nProbe SR: %.1f', sr_probe_bruce_unm), 'FontWeight', 'bold');
title('Bruce 2024');
xlim([0, t_max]); ylim([0, max_y]); set(gca, 'XTickLabel', []); box off;

% --- ROWS 2 to N: Masked PSTHs ---
for lvl_idx = 1:n_levels
    % Col 1: Original Masked
    nexttile; hold on;
    bar(t_centers, psth_orig_unm, 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    psth_m = get_psth(res_orig(lvl_idx).all_spikes);
    plot(t_centers, psth_m, 'k', 'LineWidth', 1.5);
    ylabel(sprintf('FE: %.2f%%\nMask SR: %d\nProbe SR: %.1f', ...
        target_masker_FEs(lvl_idx)*100, round(res_orig(lvl_idx).mask_sr), res_orig(lvl_idx).probe_sr));
    xlim([0, t_max]); ylim([0, max_y]); box off;
    if lvl_idx < n_levels, set(gca, 'XTickLabel', []); else, xlabel('Time (s)'); end
    
    % Col 2: Bruce Masked
    nexttile; hold on;
    bar(t_centers, psth_bruce_unm, 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    psth_m = get_psth(res_bruce(lvl_idx).all_spikes);
    plot(t_centers, psth_m, 'k', 'LineWidth', 1.5);
    ylabel(sprintf('FE: %.2f%%\nMask SR: %d\nProbe SR: %.1f', ...
        target_masker_FEs(lvl_idx)*100, round(res_bruce(lvl_idx).mask_sr), res_bruce(lvl_idx).probe_sr));
    xlim([0, t_max]); ylim([0, max_y]); box off;
    if lvl_idx < n_levels, set(gca, 'XTickLabel', []); else, xlabel('Time (s)'); end
end

disp('Simulation and plotting complete.');

%% ========================================================================
% Helper Functions
% ========================================================================

function Istim = build_safe_pulse_train(pulse, rate, duration, Fs)
    % Boundary-safe local pulse generator
    n_samples = round(duration * Fs);
    Istim = zeros(1, n_samples);
    interval = round(Fs / rate);
    onset_indices = 1 : interval : n_samples;
    
    for idx = onset_indices
        end_idx = idx + length(pulse) - 1;
        if end_idx <= n_samples
            Istim(idx:end_idx) = pulse;
        else
            Istim(idx:end) = pulse(1 : n_samples - idx + 1);
        end
    end
end

function [all_spikes, masker_sr, probe_sr] = run_model_MC(model_func, stim, Fs, NoiseAlpha, nTrials, masker_end, probe_start, masker_dur, probe_dur)
    all_spikes = [];
    masker_spikes_count = 0;
    probe_spikes_count = 0;
    
    for tr = 1:nTrials
        p_n = Library.oneonfnoise(length(stim), NoiseAlpha);
        c_n = Library.oneonfnoise(length(stim), NoiseAlpha);
        
        [~, sp, ~, ~] = model_func(stim, p_n, c_n, Fs);
        all_spikes = [all_spikes; sp(:)];
        
        masker_spikes_count = masker_spikes_count + sum(sp <= masker_end);
        probe_spikes_count  = probe_spikes_count + sum(sp >= probe_start);
    end
    
    masker_sr = masker_spikes_count / (nTrials * masker_dur);
    probe_sr  = probe_spikes_count / (nTrials * probe_dur);
end

function best_amp = calibrate_model_current(model_handle, pulse_shape, rate, calib_duration, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials)
    % Generates a short stimulus train specifically for calibration
    Istim_calib = build_safe_pulse_train(pulse_shape, rate, calib_duration, Fs);
    
    max_iters = 10; % 10 iterations of binary search
    tolerance = max(10, target_rate * 0.05); % 5% error margin or 10 sp/s min
    
    best_amp = (min_I + max_I) / 2;
    
    for iter = 1:max_iters
        test_amp = (min_I + max_I) / 2;
        input_I = Istim_calib * test_amp;
        
        total_spikes = 0;
        for tr = 1:nTrials
            p_n = Library.oneonfnoise(length(input_I), NoiseAlpha);
            c_n = Library.oneonfnoise(length(input_I), NoiseAlpha);
            [~, SpT, ~, ~, ~, ~] = model_handle(input_I, p_n, c_n, Fs);
            total_spikes = total_spikes + length(SpT);
        end
        
        current_rate = (total_spikes / nTrials) / calib_duration;
        
        if abs(current_rate - target_rate) < tolerance
            best_amp = test_amp;
            break;
        elseif current_rate < target_rate
            min_I = test_amp; % Need more current
        else
            max_I = test_amp; % Need less current
        end
        best_amp = test_amp;
    end
end