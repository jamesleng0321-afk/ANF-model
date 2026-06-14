%% ========================================================================
% Miller2011_PSTH_Biophysical.m
% Miller 2011 Fig 1 - PSTH
% Probe recovery after high-rate masker (Time-course analysis)
% Comparison: Original Model (Joshi 2017) vs Bruce (2024) Improved Model
% ========================================================================
clear; clc; close all;

%% 1. Parameters Configuration
Fs = 1e6;                       % Sampling frequency (1 MHz)
NoiseAlpha = 0.8;               % 1/f noise spectral shaping parameter
nTrials = 50;                   % Monte Carlo trials (50-100 recommended for smooth PSTH)

% Masker settings: 200 ms duration, 5000 pps
masker_duration = 0.2;          % seconds
masker_rate     = 5000;         % pps

% Probe settings: 250 ms duration, 100 pps
probe_duration  = 0.25;         % seconds
probe_rate      = 100;          % pps

% Gap between Masker offset and Probe onset
gap_duration    = 0.2e-3;       % 0.2 ms

% Stimulus levels (dB re I50)
probe_dB  = 0.0; 
masker_dB = [0.5, 0, -0.4, -0.5, -0.7]; 

% Pulse Definition: 40 us/phase symmetric biphasic, leading cathodic, 0 IPG
phase_len_us = 40;
pulse_ANF = [0, -1*ones(1, phase_len_us), 1*ones(1, phase_len_us), 0];

%% 2. Calculate Baseline Threshold (I50)
disp('Calculating Single Pulse Baseline Threshold (I50)...');
[Level_curve, Prob_curve] = Library.FindThreshold([pulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, ~] = Library.FitNeuronDynamicRange(Level_curve', Prob_curve);
I50 = muSingle;
fprintf('I50 found: %.4f uA\n', I50 * 1e6);

%% 3. Construct Base Stimuli Vectors
% [核心修复]: 使用本地的安全脉冲生成器，绕过 Experiment.stim_PulseTrain 的数组溢出 Bug
masker_base = build_safe_pulse_train(pulse_ANF, masker_rate, masker_duration, Fs);
probe_base  = build_safe_pulse_train(pulse_ANF, probe_rate, probe_duration, Fs);
gap_zeros   = zeros(1, round(gap_duration * Fs));

% Define timing boundaries for analysis
masker_end  = length(masker_base) / Fs;
probe_start = masker_end + gap_duration;
t_max       = probe_start + probe_duration;

%% 4. Unmasked Baselines (Probe Only)
fprintf('\n--- Simulating Unmasked Probe Baselines ---\n');
amp_probe = I50 * 10^(probe_dB/20);
stim_unmasked = [zeros(1, length(masker_base) + length(gap_zeros)), probe_base * amp_probe];

[spikes_orig_unmasked, ~, sr_probe_orig_unm] = run_model_MC(@Model_PulseTrain, stim_unmasked, Fs, NoiseAlpha, nTrials, masker_end, probe_start, masker_duration, probe_duration);
[spikes_bruce_unmasked, ~, sr_probe_bruce_unm] = run_model_MC(@Model_PulseTrain_Bruce, stim_unmasked, Fs, NoiseAlpha, nTrials, masker_end, probe_start, masker_duration, probe_duration);

%% 5. Masked Simulation Loop
fprintf('\n--- Running Masked Accomodation Simulations ---\n');
n_levels = length(masker_dB);

% Data structures to store results
res_orig  = struct('all_spikes', {}, 'mask_sr', {}, 'probe_sr', {});
res_bruce = struct('all_spikes', {}, 'mask_sr', {}, 'probe_sr', {});

for lvl_idx = 1:n_levels
    level = masker_dB(lvl_idx);
    fprintf('Processing Masker Level: %.1f dB re I50...\n', level);
    
    amp_mask = I50 * 10^(level/20);
    stim_masked = [masker_base * amp_mask, gap_zeros, probe_base * amp_probe];
    
    % Run Original Model
    [sp, m_sr, p_sr] = run_model_MC(@Model_PulseTrain, stim_masked, Fs, NoiseAlpha, nTrials, masker_end, probe_start, masker_duration, probe_duration);
    res_orig(lvl_idx).all_spikes = sp; res_orig(lvl_idx).mask_sr = m_sr; res_orig(lvl_idx).probe_sr = p_sr;
    
    % Run Bruce Model
    [sp, m_sr, p_sr] = run_model_MC(@Model_PulseTrain_Bruce, stim_masked, Fs, NoiseAlpha, nTrials, masker_end, probe_start, masker_duration, probe_duration);
    res_bruce(lvl_idx).all_spikes = sp; res_bruce(lvl_idx).mask_sr = m_sr; res_bruce(lvl_idx).probe_sr = p_sr;
end

%% 6. Plotting PSTHs
% Configuration for PSTH binning
bin_width = 2e-3; % 2 ms bins for PSTH calculation
edges = 0 : bin_width : t_max;
t_centers = edges(1:end-1) + bin_width/2;

fig = figure("Name", "Miller2011_PSTH_Biophysical_Comparison", "DefaultAxesFontSize", 11, "Color", "w");
fig.Position(3:4) = [1000, 900];
% 2 Columns: Original (Left) vs Bruce (Right)
tiledlayout(n_levels + 1, 2, "TileSpacing", "compact", "Padding", "compact");

% Helper function for PSTH binning
get_psth = @(spikes) histcounts(spikes, edges) / (nTrials * bin_width);

% --- ROW 1: Unmasked Baselines ---
psth_orig_unm = get_psth(spikes_orig_unmasked);
psth_bruce_unm = get_psth(spikes_bruce_unmasked);
max_y = max([max(psth_orig_unm), max(psth_bruce_unm)]) * 1.2;

% Col 1: Original Unmasked
nexttile; hold on;
bar(t_centers, psth_orig_unm, 1, 'FaceColor', '#0072BD', 'EdgeColor', 'none');
ylabel(sprintf('Unmasked\nProbe SR: %.1f', sr_probe_orig_unm), 'FontWeight', 'bold');
title('Original Model (Joshi 2017)');
xlim([0, t_max]); ylim([0, max_y]); set(gca, 'XTickLabel', []); box off;

% Col 2: Bruce Unmasked
nexttile; hold on;
bar(t_centers, psth_bruce_unm, 1, 'FaceColor', '#D95319', 'EdgeColor', 'none');
ylabel(sprintf('Unmasked\nProbe SR: %.1f', sr_probe_bruce_unm), 'FontWeight', 'bold');
title('Bruce (2024) Improved Model');
xlim([0, t_max]); ylim([0, max_y]); set(gca, 'XTickLabel', []); box off;

% --- ROWS 2 to N: Masked PSTHs ---
for lvl_idx = 1:n_levels
    % Col 1: Original Masked
    nexttile; hold on;
    bar(t_centers, psth_orig_unm, 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    psth_m = get_psth(res_orig(lvl_idx).all_spikes);
    plot(t_centers, psth_m, 'k', 'LineWidth', 1.5);
    ylabel(sprintf('Mask: %.1f dB\nMask SR: %d\nProbe SR: %.1f', ...
        masker_dB(lvl_idx), round(res_orig(lvl_idx).mask_sr), res_orig(lvl_idx).probe_sr));
    xlim([0, t_max]); ylim([0, max_y]); box off;
    if lvl_idx < n_levels, set(gca, 'XTickLabel', []); else, xlabel('Time [s]'); end
    
    % Col 2: Bruce Masked
    nexttile; hold on;
    bar(t_centers, psth_bruce_unm, 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    psth_m = get_psth(res_bruce(lvl_idx).all_spikes);
    plot(t_centers, psth_m, 'k', 'LineWidth', 1.5);
    ylabel(sprintf('Mask: %.1f dB\nMask SR: %d\nProbe SR: %.1f', ...
        masker_dB(lvl_idx), round(res_bruce(lvl_idx).mask_sr), res_bruce(lvl_idx).probe_sr));
    xlim([0, t_max]); ylim([0, max_y]); box off;
    if lvl_idx < n_levels, set(gca, 'XTickLabel', []); else, xlabel('Time [s]'); end
end

%% ========================================================================
% Helper Functions
% ========================================================================

function Istim = build_safe_pulse_train(pulse, rate, duration, Fs)
    % 边界安全的局部脉冲生成器
    n_samples = round(duration * Fs);
    Istim = zeros(1, n_samples);
    interval = round(Fs / rate);
    onset_indices = 1 : interval : n_samples;
    
    for idx = onset_indices
        end_idx = idx + length(pulse) - 1;
        if end_idx <= n_samples
            Istim(idx:end_idx) = pulse;
        else
            % 如果脉冲在时间窗边缘溢出，则强制截断以保护数组维度
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