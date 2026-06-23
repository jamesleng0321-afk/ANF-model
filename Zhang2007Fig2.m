%% ========================================================================
% Evaluate_Adaptation_BouletStyle.m
% Spike Rate Adaptation evaluation for Joshi 2017 and Bruce 2024 models.
% Replicates the style of Boulet 2017 Fig 7 using REAL empirical data
% from Zhang 2007 Fig 2.
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Load Empirical Data (Zhang 2007)
% ------------------------------------------------------------------------
fprintf('Loading empirical data from Zhang_2007_fig2.csv...\n');
try
    zhang_data = readtable('Zhang_2007_fig2.csv');
catch
    error('Could not find Zhang_2007_fig2.csv. Please ensure it is in the current directory.');
end

% Extract time points (ms)
exp_time = zhang_data.x;

% Define the column names for the 5000 pps experimental data
% Note: Using the exact column names from the CSV (including the typo 'adapation' in the CSV)
data_columns = {'adaptation_5000pps_0_35mA', ...
                'adapation_5000pps_0_38mA', ...
                'adaptation_5000pps_0_5mA'};

% Extract the exact onset spike rates (the first data point of each column)
% to use as calibration targets so models and data start at the same point.
target_rates = [zhang_data.(data_columns{1})(1), ...
                zhang_data.(data_columns{2})(1), ...
                zhang_data.(data_columns{3})(1)];

%% ------------------------------------------------------------------------
% 2. Global Parameters
% ------------------------------------------------------------------------
stim_rate = 5000;           % Stimulus rate in pulses per second (pps)
duration_s = 300e-3;        % Full simulation duration: 300 ms
onset_window = 12e-3;       % Onset analysis window: 12 ms
bin_width = 12e-3;          % Bin size for spike rate calculation: 12 ms

nTrials_Calib = 20;         % Number of trials for the calibration phase
nTrials_Sim = 50;           % Number of trials for the full simulation phase (increase for smoother curves)

Fs = 1e6;
NoiseAlpha = 0.8;
% Biphasic pulse: 40us negative, 0us IPG, 40us positive
SinglePulse = [0, -1*ones(1,40), 1*ones(1,40), 0]; 

%% ------------------------------------------------------------------------
% 3. Auto-Calibration Phase (Binary Search for empirical Onset Rates)
% ------------------------------------------------------------------------
fprintf('\n=== Phase 1: Auto-Calibrating Currents to Match Empirical Onsets ===\n');

calibrated_amps_Joshi = zeros(1, length(target_rates));
calibrated_amps_Bruce = zeros(1, length(target_rates));

% Define search boundaries for current (Amperes)
min_I = 0.5e-3; 
max_I = 2.5e-3; 

for i = 1:length(target_rates)
    target_rate = target_rates(i);
    target_FE = target_rate / stim_rate; 
    fprintf('\nTarget Onset Rate: %.1f spikes/s (FE: %.1f%%)\n', target_rate, target_FE*100);
    
    % Calibrate Joshi 2017
    calibrated_amps_Joshi(i) = calibrate_model_current(@Model_PulseTrain, ...
        SinglePulse, stim_rate, onset_window, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials_Calib);
    fprintf('  -> Joshi 2017 calibrated current: %.3f mA\n', calibrated_amps_Joshi(i)*1000);
    
    % Calibrate Bruce 2024
    calibrated_amps_Bruce(i) = calibrate_model_current(@Model_PulseTrain_Bruce, ...
        SinglePulse, stim_rate, onset_window, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials_Calib);
    fprintf('  -> Bruce 2024 calibrated current: %.3f mA\n', calibrated_amps_Bruce(i)*1000);
end

%% ------------------------------------------------------------------------
% 4. Full Simulation Phase (300 ms)
% ------------------------------------------------------------------------
fprintf('\n=== Phase 2: Running Full 300 ms Simulations ===\n');

% Define time bins for the adaptation curve
time_edges = 0 : bin_width : duration_s;
time_centers = time_edges(1:end-1) + bin_width/2;

rate_Joshi_all = cell(1, length(target_rates));
rate_Bruce_all = cell(1, length(target_rates));

for i = 1:length(target_rates)
    fprintf('Running simulations for target %.1f spikes/s...\n', target_rates(i));
    
    % Generate full 300ms stimulus train
    Istim_full = Experiment.stim_PulseTrain(SinglePulse, stim_rate, 100, 0, duration_s, Fs);
    
    input_Joshi = Istim_full * calibrated_amps_Joshi(i);
    input_Bruce = Istim_full * calibrated_amps_Bruce(i);
    
    all_spikes_Joshi = [];
    all_spikes_Bruce = [];
    
    for tr = 1:nTrials_Sim
        p_noise = Library.oneonfnoise(length(input_Joshi), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(input_Joshi), NoiseAlpha);
        
        % Run Joshi
        [~, SpTimes_J, ~, ~, ~, ~] = Model_PulseTrain(input_Joshi, p_noise, c_noise, Fs);
        all_spikes_Joshi = [all_spikes_Joshi; SpTimes_J(:)];
        
        % Run Bruce
        [~, SpTimes_B, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_Bruce, p_noise, c_noise, Fs);
        all_spikes_Bruce = [all_spikes_Bruce; SpTimes_B(:)];
    end
    
    % Calculate spike rate in each 12ms bin (spikes / sec)
    counts_J = histcounts(all_spikes_Joshi, time_edges);
    rate_Joshi_all{i} = (counts_J / nTrials_Sim) / bin_width;
    
    counts_B = histcounts(all_spikes_Bruce, time_edges);
    rate_Bruce_all{i} = (counts_B / nTrials_Sim) / bin_width;
end

%% ------------------------------------------------------------------------
% 5. Plotting (Boulet 2017 Fig 7 Style)
% ------------------------------------------------------------------------
fprintf('\n=== Phase 3: Plotting Results ===\n');

fig = figure("Name", "Spike Rate Adaptation", "DefaultAxesFontSize", 12);
fig.Position(3:4) = [1000 800];
tl = tiledlayout(length(target_rates), 1, "TileSpacing", "compact", "Padding", "compact");
title(tl, sprintf('Spike Rate Adaptation (%d pps)', stim_rate), 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:length(target_rates)
    nexttile; hold on;
    
    % A. Plot Empirical Data (Zhang 2007) as scatter points
    exp_rate = zhang_data.(data_columns{i});
    l_exp = plot(exp_time, exp_rate, 'ko', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 6);
    
    % B. Plot Model Predictions as solid lines
    l_joshi = plot(time_centers*1000, rate_Joshi_all{i}, '-', 'LineWidth', 2, 'Color', "#D95319");
    l_bruce = plot(time_centers*1000, rate_Bruce_all{i}, '-', 'LineWidth', 2, 'Color', "#0072BD");
    
    % Formatting
    xlim([0 duration_s*1000]);
    % Dynamically set Y limit based on the maximum onset rate for better visibility
    ylim([0 max(100, target_rates(i) * 1.2)]); 
    grid on; box off;
    
    % Construct Y Label identifying the specific target onset
    ylabel(sprintf('Spike Rate (spikes/s)\nOnset ~%.0f sp/s', target_rates(i)), 'FontWeight', 'bold');
    
    if i == length(target_rates)
        xlabel('Time (ms)', 'FontWeight', 'bold');
    end
    
    % Add legend only to the first subplot
    if i == 1
        lg = legend([l_exp, l_joshi, l_bruce], ["Zhang 2007 (Cat Data)", "Joshi 2017", "Bruce 2024"]);
        lg.Orientation = 'horizontal';
        lg.Location = 'northeast';
    end
end

disp('Simulation and plotting complete.');


%% ========================================================================
% Helper Function: Binary Search Calibration for Target Onset Rate
% ========================================================================
function best_amp = calibrate_model_current(model_handle, pulse_shape, rate, calib_duration, target_rate, min_I, max_I, Fs, NoiseAlpha, nTrials)
    % Generates a short stimulus train specifically for calibration (e.g., 12 ms)
    Istim_calib = Experiment.stim_PulseTrain(pulse_shape, rate, 100, 0, calib_duration, Fs);
    
    max_iters = 10; % 10 iterations of binary search
    tolerance = target_rate * 0.05; % Allow 5% error margin
    
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
            min_I = test_amp; % Need more current to increase spike rate
        else
            max_I = test_amp; % Need less current to decrease spike rate
        end
        best_amp = test_amp;
    end
end