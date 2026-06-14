%% ========================================================================
% Miller2011_Accommodation_Biophysical_Models.m
% Miller 2011 Fig 2 - Accommodation (Probe recovery after high-rate masker)
% Comparison: Original Model (Joshi 2017) vs Bruce (2024) Improved Model
% ========================================================================
clear; clc; close all;

%% 1. Load Experimental Data
load Miller_2011_fig2.mat % Ensure this file is in your path

%% 2. Parameters Configuration
Fs = 1e6;                       % 1 MHz sampling rate
NoiseAlpha = 0.8;               % 1/f noise shaping

target_probe_FE   = 0.7;        % Probe level fixed to achieve ~70% FE
threshold_mask_FE = 0.005;      % Masker threshold reference (~0.5% FE)

masker_rates    = [250 5000];   % pps
probe_rate      = 100;          % pps
probe_duration  = 0.3;          % seconds
masker_duration = 0.3;          % seconds

levels_dB = -4:1:5;             % Masker levels (relative to threshold)

% Pulse Definition: 40 us/phase symmetric biphasic, leading cathodic, 0 IPG
phase_len_us = 40;
pulse_ANF = [0, -1*ones(1, phase_len_us), 1*ones(1, phase_len_us), 0];

%% 3. Find Baseline Thresholds & Amplitudes
disp('Calculating single pulse dynamic range to estimate levels...');
[Level_curve, Prob_curve] = Library.FindThreshold([pulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, sigmaSingle]   = Library.FitNeuronDynamicRange(Level_curve', Prob_curve);

% Original ANF Amplitude Mapping
get_level_ANF = @(p) muSingle + sigmaSingle * sqrt(2) * erfinv(2*p - 1);
probe_amp_ANF = get_level_ANF(target_probe_FE);
mask_amp_ANF  = get_level_ANF(threshold_mask_FE);

% Bruce 2024 Amplitude Mapping (Sigma scaled by 1/3)
get_level_Bruce = @(p) muSingle + (sigmaSingle / 3) * sqrt(2) * erfinv(2*p - 1);
probe_amp_Bruce = get_level_Bruce(target_probe_FE);
mask_amp_Bruce  = get_level_Bruce(threshold_mask_FE);

%% 4. Stimulus Construction Setup
% Construct normalized base pulse trains
mask_250_base  = Experiment.stim_PulseTrain(pulse_ANF, masker_rates(1), 100, 0, masker_duration, Fs);
mask_5000_base = Experiment.stim_PulseTrain(pulse_ANF, masker_rates(2), 100, 0, masker_duration, Fs);
probe_base     = Experiment.stim_PulseTrain(pulse_ANF, probe_rate, 100, 0, probe_duration, Fs);

num_probe_pulses = floor(probe_duration * probe_rate);
nTrials = 20; % Monte Carlo trials (Increase to 50 for smoother curves)

%% 5. Unmasked Baselines
fprintf('Calculating unmasked probe responses...\n');
unmasked_FE_ANF   = compute_probe_FE(@Model_PulseTrain, [zeros(1, length(mask_250_base)), probe_base * probe_amp_ANF], masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);
unmasked_FE_Bruce = compute_probe_FE(@Model_PulseTrain_Bruce, [zeros(1, length(mask_250_base)), probe_base * probe_amp_Bruce], masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);

%% 6. Main Monte Carlo Simulation Loop
n_levels = length(levels_dB);
recov250_ANF  = zeros(n_levels,1); recov5000_ANF  = zeros(n_levels,1);
recov250_Bruce= zeros(n_levels,1); recov5000_Bruce= zeros(n_levels,1);

fprintf('Running accommodation simulation...\n');
for lvl_idx = 1:n_levels
    fac = 10^(levels_dB(lvl_idx)/20);
    fprintf('  Processing Level: %d dB...\n', levels_dB(lvl_idx));
    
    % --- 250 pps Masker ---
    stim_250_ANF   = [mask_250_base * mask_amp_ANF * fac, probe_base * probe_amp_ANF];
    stim_250_Bruce = [mask_250_base * mask_amp_Bruce * fac, probe_base * probe_amp_Bruce];
    
    recov250_ANF(lvl_idx)   = compute_probe_FE(@Model_PulseTrain, stim_250_ANF, masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);
    recov250_Bruce(lvl_idx) = compute_probe_FE(@Model_PulseTrain_Bruce, stim_250_Bruce, masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);

    % --- 5000 pps Masker ---
    stim_5000_ANF   = [mask_5000_base * mask_amp_ANF * fac, probe_base * probe_amp_ANF];
    stim_5000_Bruce = [mask_5000_base * mask_amp_Bruce * fac, probe_base * probe_amp_Bruce];
    
    recov5000_ANF(lvl_idx)   = compute_probe_FE(@Model_PulseTrain, stim_5000_ANF, masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);
    recov5000_Bruce(lvl_idx) = compute_probe_FE(@Model_PulseTrain_Bruce, stim_5000_Bruce, masker_duration, num_probe_pulses, Fs, NoiseAlpha, nTrials);
end

%% 7. Normalize by Unmasked Response (Recovery Ratio)
recov250_ANF  = recov250_ANF  / max(unmasked_FE_ANF, 1e-4);
recov5000_ANF = recov5000_ANF / max(unmasked_FE_ANF, 1e-4);

recov250_Bruce  = recov250_Bruce  / max(unmasked_FE_Bruce, 1e-4);
recov5000_Bruce = recov5000_Bruce / max(unmasked_FE_Bruce, 1e-4);

%% 8. Plotting Results
fig = figure("Name","Miller2011_Accommodation_Comparison", "DefaultAxesFontSize",13, "Color","w");
fig.Position(3:4) = [1100 450];
tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

% ---- 250 pps Subplot ----
nexttile; hold on;
scatter(accommodation_Miller2011.Masker250.Raw.MaskerLevel, ...
        accommodation_Miller2011.Masker250.Raw.Recovery, ...
        'o','MarkerEdgeColor','#808080','MarkerEdgeAlpha',0.4, 'DisplayName', 'Raw Data');
plot(accommodation_Miller2011.Masker250.Medians.MaskerLevel, ...
     accommodation_Miller2011.Masker250.Medians.Recovery, ...
     "kx-","LineWidth",2,"MarkerSize",8, 'DisplayName', 'Median Data');
plot(levels_dB, recov250_ANF, "d-", "Color", "#0072BD", "LineWidth",2, 'DisplayName', 'Original Model');
plot(levels_dB, recov250_Bruce, "^-", "Color", "#D95319", "LineWidth",2, 'DisplayName', 'Bruce (2024) Model');

xlabel("Masker level [dB re thr]"); ylabel("Probe recovery ratio");
title("250 pps Masker"); grid on; box off;

% ---- 5000 pps Subplot ----
nexttile; hold on;
scatter(accommodation_Miller2011.Masker5000.Raw.MaskerLevel, ...
        accommodation_Miller2011.Masker5000.Raw.Recovery, ...
        'o','MarkerEdgeColor','#808080','MarkerEdgeAlpha',0.4, 'DisplayName', 'Raw Data');
plot(accommodation_Miller2011.Masker5000.Medians.MaskerLevel, ...
     accommodation_Miller2011.Masker5000.Medians.Recovery, ...
     "kx-","LineWidth",2,"MarkerSize",8, 'DisplayName', 'Median Data');
plot(levels_dB, recov5000_ANF, "d-", "Color", "#0072BD", "LineWidth",2, 'DisplayName', 'Original Model');
plot(levels_dB, recov5000_Bruce, "^-", "Color", "#D95319", "LineWidth",2, 'DisplayName', 'Bruce (2024) Model');

xlabel("Masker level [dB re thr]"); ylabel("Probe recovery ratio");
title("5000 pps Masker"); grid on; box off;

lg = legend("Orientation","horizontal");
lg.Layout.Tile = "north";

%% ========================================================================
% Helper Functions
% ========================================================================

function probe_FE = compute_probe_FE(model_func, stim_input, probe_start_time, num_probe_pulses, Fs, NoiseAlpha, nTrials)
    total_probe_spikes = 0;
    
    for tr = 1:nTrials
        p_noise = Library.oneonfnoise(length(stim_input), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(stim_input), NoiseAlpha);
        
        [~, SpTimes, ~, ~] = model_func(stim_input, p_noise, c_noise, Fs);
        
        % Only count spikes that occur after the masker duration
        probe_spikes_in_trial = sum(SpTimes > probe_start_time);
        total_probe_spikes = total_probe_spikes + probe_spikes_in_trial;
    end
    
    probe_FE = total_probe_spikes / (nTrials * num_probe_pulses);
end