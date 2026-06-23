%% ========================================================================
% Evaluate_Dynes1996_Facilitation.m
% Reproduction of Dynes 1996 Fig 4-1: Facilitation and Masking
% Probe threshold relative to single-pulse threshold as a function of MPI.
% Comparison: Joshi 2017 vs Bruce 2024
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Load Experimental Data (Dynes 1996)
% ------------------------------------------------------------------------
try
    measured_data = readtable("Dynes_1996_fig4-1.csv");
    has_data = true;
    fprintf('Loaded Dynes_1996_fig4-1.csv successfully.\n');
catch
    warning('Dynes_1996_fig4-1.csv not found. Using default MPI values for simulation.');
    has_data = false;
end

masker_counts = [1, 4, 40];      % Number of maskers used in the paper
phase_len_us  = 100;             % 100 us phase duration (matches Dynes/aLIFP params)

% Fs and Noise
Fs = 1e6;
NoiseAlpha = 0.8;
nTrials = 20; % Monte Carlo trials for threshold estimation

% Define 100us biphasic pulse (cathodic first)
SinglePulse = [0, -1*ones(1, phase_len_us), 1*ones(1, phase_len_us), 0];

%% ------------------------------------------------------------------------
% 2. Compute Baseline Single-Pulse Threshold (I50)
% ------------------------------------------------------------------------
fprintf('=== Phase 1: Calculating Baseline Single-Pulse Thresholds ===\n');

% We use a targeted sweep to find the exact I50 for both models
amp_sweep_baseline = linspace(0.4e-3, 1.5e-3, 15); % 0.4 mA to 1.5 mA

I50_single_Joshi = get_probe_threshold(@Model_PulseTrain, 0, 0, 0, SinglePulse, Fs, NoiseAlpha, nTrials, amp_sweep_baseline);
I50_single_Bruce = get_probe_threshold(@Model_PulseTrain_Bruce, 0, 0, 0, SinglePulse, Fs, NoiseAlpha, nTrials, amp_sweep_baseline);

fprintf('  -> Joshi 2017 Baseline I50: %.3f mA\n', I50_single_Joshi * 1000);
fprintf('  -> Bruce 2024 Baseline I50: %.3f mA\n', I50_single_Bruce * 1000);

%% ------------------------------------------------------------------------
% 3. Masker Amplitude (89% of Single-Pulse Threshold)
% ------------------------------------------------------------------------
% The experiment uses sub-threshold maskers fixed at 1 dB below threshold
masker_amp_Joshi = 0.89 * I50_single_Joshi;
masker_amp_Bruce = 0.89 * I50_single_Bruce;

%% ------------------------------------------------------------------------
% 4. Main Simulation Loop
% ------------------------------------------------------------------------
fprintf('\n=== Phase 2: Running Masker-Probe Simulations ===\n');

% Set up figure (1x3 tiled layout for 1, 4, and 40 maskers)
fig = figure("Name", "Dynes1996_Facilitation_Masking", "DefaultAxesFontSize", 12);
fig.Position(3:4) = [1200 450];
tl = tiledlayout(1, length(masker_counts), "TileSpacing", "compact", "Padding", "compact");
title(tl, 'Probe Threshold Shift vs Masker-Probe Interval (Dynes 1996)', 'FontSize', 16, 'FontWeight', 'bold');

for m_ind = 1:length(masker_counts)
    num_masker = masker_counts(m_ind);
    fprintf('\nSimulating condition: %d Masker(s)...\n', num_masker);
    
    % Determine MPIs (Masker-Probe Intervals) to test
    if has_data
        mpi_data = measured_data.(sprintf('x%d', num_masker)) * 1e-3; % ms to s
        ratio_db = measured_data.(sprintf('y%d', num_masker));
        
        valid_idx = ~isnan(mpi_data);
        mpis_s = mpi_data(valid_idx);
        data_ratio_linear = 10.^(ratio_db(valid_idx) / 20); % Convert dB to linear ratio
    else
        % Default test points if CSV is missing
        mpis_s = [0.5, 0.8, 1, 1.5, 2, 4, 8, 16] * 1e-3; 
    end
    
    I50_ratio_Joshi = zeros(1, length(mpis_s));
    I50_ratio_Bruce = zeros(1, length(mpis_s));
    
    % Sweep range for probe threshold relative to baseline (from 40% to 150%)
    amp_multipliers = linspace(0.4, 1.5, 15);
    sweep_Joshi = amp_multipliers * I50_single_Joshi;
    sweep_Bruce = amp_multipliers * I50_single_Bruce;
    
    for i = 1:length(mpis_s)
        current_mpi = mpis_s(i);
        
        % Evaluate Joshi 2017
        probe_I50_J = get_probe_threshold(@Model_PulseTrain, num_masker, masker_amp_Joshi, current_mpi, SinglePulse, Fs, NoiseAlpha, nTrials, sweep_Joshi);
        I50_ratio_Joshi(i) = probe_I50_J / I50_single_Joshi;
        
        % Evaluate Bruce 2024
        probe_I50_B = get_probe_threshold(@Model_PulseTrain_Bruce, num_masker, masker_amp_Bruce, current_mpi, SinglePulse, Fs, NoiseAlpha, nTrials, sweep_Bruce);
        I50_ratio_Bruce(i) = probe_I50_B / I50_single_Bruce;
    end
    
    % ---- Plotting ----
    nexttile; hold on;
    
    % Plot 1.0 reference line (No shift)
    plot([0, max(mpis_s)*1e3], [1, 1], 'k--', 'LineWidth', 1);
    
    % Plot Experimental Data
    if has_data
        l_exp = plot(mpis_s*1e3, data_ratio_linear, 'x', 'MarkerSize', 8, 'LineWidth', 2, 'Color', '#7E2F8E');
    end
    
    % Plot Joshi 2017
    l_joshi = plot(mpis_s*1e3, I50_ratio_Joshi, '-o', 'LineWidth', 2, 'Color', "#D95319", 'MarkerFaceColor', "#D95319");
    
    % Plot Bruce 2024
    l_bruce = plot(mpis_s*1e3, I50_ratio_Bruce, '-s', 'LineWidth', 2, 'Color', "#0072BD", 'MarkerFaceColor', "#0072BD");
    
    % Formatting
    grid on; box off;
    ylim([0.4 1.4]);
    xlim([0 max(mpis_s)*1e3 + 1]);
    title(sprintf('%d Masker(s)', num_masker));
    
    if m_ind == 1
        ylabel('Threshold Ratio (I_{50, probe} / I_{50, single})', 'FontWeight', 'bold');
    end
    if m_ind == 2
        xlabel('Masker-Probe Interval [ms]', 'FontWeight', 'bold');
    end
    
    % Legend on first tile
    if m_ind == 1
        if has_data
            legend([l_exp, l_joshi, l_bruce], {"Dynes 1996", "Joshi 2017", "Bruce 2024"}, 'Location', 'southeast');
        else
            legend([l_joshi, l_bruce], {"Joshi 2017", "Bruce 2024"}, 'Location', 'southeast');
        end
    end
end

disp('Simulation complete.');


%% ========================================================================
% Helper Functions
% ========================================================================

function [stim, probe_onset_s] = build_dynes_stimulus(num_maskers, masker_amp, mpi_s, probe_amp, pulse, Fs)
    % Builds the stimulus vector. Maskers are spaced by 1 ms (1000 pps).
    masker_ipi = 1e-3; 
    
    if num_maskers > 0
        last_masker_onset = (num_maskers - 1) * masker_ipi;
    else
        last_masker_onset = 0;
        mpi_s = 0; % Force 0 MPI if it's just a single baseline pulse
    end
    
    probe_onset_s = last_masker_onset + mpi_s;
    total_duration = probe_onset_s + 5e-3; % Add 5 ms padding after probe to capture spikes
    
    stim = zeros(1, round(total_duration * Fs));
    pulse_len = length(pulse);
    
    % Insert maskers
    for k = 1:num_maskers
        idx = round((k-1) * masker_ipi * Fs) + 1;
        stim(idx : idx + pulse_len - 1) = pulse * masker_amp;
    end
    
    % Insert probe
    p_idx = round(probe_onset_s * Fs) + 1;
    stim(p_idx : p_idx + pulse_len - 1) = pulse * probe_amp;
end

function threshold = get_probe_threshold(model_handle, num_maskers, masker_amp, mpi_s, pulse, Fs, NoiseAlpha, nTrials, test_amps)
    % Sweeps through test amplitudes to find the 50% firing probability of the probe
    probs = zeros(1, length(test_amps));
    
    for a_idx = 1:length(test_amps)
        [stim, probe_start_s] = build_dynes_stimulus(num_maskers, masker_amp, mpi_s, test_amps(a_idx), pulse, Fs);
        
        spikes_count = 0;
        for tr = 1:nTrials
            p_n = Library.oneonfnoise(length(stim), NoiseAlpha);
            c_n = Library.oneonfnoise(length(stim), NoiseAlpha);
            [~, SpT, ~, ~, ~, ~] = model_handle(stim, p_n, c_n, Fs);
            
            % Count trial as "spiked" if at least one spike occurs in the 4ms window following the probe
            if any(SpT >= probe_start_s & SpT <= probe_start_s + 4e-3)
                spikes_count = spikes_count + 1;
            end
        end
        probs(a_idx) = spikes_count / nTrials;
    end
    
    % Extract the 50% threshold using robust linear interpolation
    threshold = extract_50_threshold(test_amps, probs);
end

function threshold = extract_50_threshold(amps, probs)
    % Smooth probabilities slightly to handle Monte Carlo stochastic noise
    smoothed_probs = movmean(probs, 3);
    
    idx = find(smoothed_probs >= 0.5, 1);
    
    if isempty(idx)
        threshold = amps(end); % Maxed out, return highest tested amplitude
    elseif idx == 1
        threshold = amps(1);   % Very sensitive, return lowest tested amplitude
    else
        % Linear interpolation between the two points crossing 0.5
        x1 = smoothed_probs(idx-1); x2 = smoothed_probs(idx);
        y1 = amps(idx-1); y2 = amps(idx);
        if x2 == x1
            threshold = y1;
        else
            threshold = y1 + (0.5 - x1) * (y2 - y1) / (x2 - x1);
        end
    end
end