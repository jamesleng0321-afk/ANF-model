%% ========================================================================
% Evaluate_Decrement_MultiRate.m
% Evaluates Spike Rate Decrement vs Onset Rate across multiple pulse rates.
% Replicates the full multi-panel style of Zhang 2007 Fig 5 / Boulet 17 Fig 8.
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Global Parameters
% ------------------------------------------------------------------------
% The fundamental stimulus rates used in the original studies
pulse_rates = [250, 1000, 5000]; 

duration_s = 300e-3;        % Full simulation duration: 300 ms
onset_window = [0, 12e-3];  % Onset analysis window: 0 to 12 ms
steady_window = [200e-3, 300e-3]; % Steady-state window: 200 to 300 ms

nTrials = 20;               % Number of Monte Carlo trials per current level

Fs = 1e6;
NoiseAlpha = 0.8;
% Biphasic pulse: 40us negative, 0us IPG, 40us positive
SinglePulse = [0, -1*ones(1,40), 1*ones(1,40), 0]; 

% Define a range of current amplitudes to sweep (Unit: mA)
% A broader sweep to ensure we hit high onset rates across all frequencies
amps_mA = linspace(0.5, 2.5, 12); 

%% ------------------------------------------------------------------------
% 2. Setup Figure and Tiled Layout
% ------------------------------------------------------------------------
fig = figure("Name", "Spike Rate Decrement - Multi Rate", "DefaultAxesFontSize", 12);
% Make the figure wide to accommodate 3 subplots side-by-side
fig.Position(3:4) = [1200 450]; 
tl = tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");
title(tl, 'Spike Rate Decrement vs Onset Rate', 'FontSize', 16, 'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% 3. Main Loop: Iterate over Pulse Rates
% ------------------------------------------------------------------------
for r = 1:length(pulse_rates)
    stim_rate = pulse_rates(r);
    fprintf('\n=== Simulating Spike Rate Decrement for %d pps ===\n', stim_rate);
    
    % Pre-allocate arrays for current rate
    onset_rates_Joshi = zeros(1, length(amps_mA));
    steady_rates_Joshi = zeros(1, length(amps_mA));
    
    onset_rates_Bruce = zeros(1, length(amps_mA));
    steady_rates_Bruce = zeros(1, length(amps_mA));
    
    % Generate the 300ms stimulus train structure for the current rate
    Istim_base = Experiment.stim_PulseTrain(SinglePulse, stim_rate, 100, 0, duration_s, Fs);
    
    % ---- Sweep Amplitudes ----
    for i = 1:length(amps_mA)
        current_amp = amps_mA(i) * 1e-3; % Convert to Amperes
        input_I = Istim_base * current_amp;
        
        % Accumulators for Joshi
        spikes_onset_J = 0;
        spikes_steady_J = 0;
        
        % Accumulators for Bruce
        spikes_onset_B = 0;
        spikes_steady_B = 0;
        
        for tr = 1:nTrials
            % Synchronized noise for both models
            p_noise = Library.oneonfnoise(length(input_I), NoiseAlpha);
            c_noise = Library.oneonfnoise(length(input_I), NoiseAlpha);
            
            % Run Joshi
            [~, SpT_J, ~, ~, ~, ~] = Model_PulseTrain(input_I, p_noise, c_noise, Fs);
            spikes_onset_J = spikes_onset_J + sum(SpT_J >= onset_window(1) & SpT_J <= onset_window(2));
            spikes_steady_J = spikes_steady_J + sum(SpT_J >= steady_window(1) & SpT_J <= steady_window(2));
            
            % Run Bruce
            [~, SpT_B, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_I, p_noise, c_noise, Fs);
            spikes_onset_B = spikes_onset_B + sum(SpT_B >= onset_window(1) & SpT_B <= onset_window(2));
            spikes_steady_B = spikes_steady_B + sum(SpT_B >= steady_window(1) & SpT_B <= steady_window(2));
        end
        
        % Calculate mean rates in spikes/s
        onset_rates_Joshi(i) = (spikes_onset_J / nTrials) / (onset_window(2) - onset_window(1));
        steady_rates_Joshi(i) = (spikes_steady_J / nTrials) / (steady_window(2) - steady_window(1));
        
        onset_rates_Bruce(i) = (spikes_onset_B / nTrials) / (onset_window(2) - onset_window(1));
        steady_rates_Bruce(i) = (spikes_steady_B / nTrials) / (steady_window(2) - steady_window(1));
    end
    
    % Calculate Decrement: Onset Rate - Steady-State Rate
    decrement_Joshi = onset_rates_Joshi - steady_rates_Joshi;
    decrement_Bruce = onset_rates_Bruce - steady_rates_Bruce;
    
    % ---- Plotting the current rate panel ----
    nexttile; hold on;
    
    % Plot the y = x reference line (representing 100% adaptation)
    max_val = max([onset_rates_Joshi, onset_rates_Bruce, 100]); % Ensure at least 100 scale
    l_ref = plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.5);
    
    % Plot Joshi 2017
    l_joshi = plot(onset_rates_Joshi, decrement_Joshi, '-o', 'LineWidth', 2, ...
        'MarkerSize', 5, 'Color', "#D95319", 'MarkerFaceColor', "#D95319");
    
    % Plot Bruce 2024
    l_bruce = plot(onset_rates_Bruce, decrement_Bruce, '-s', 'LineWidth', 2, ...
        'MarkerSize', 5, 'Color', "#0072BD", 'MarkerFaceColor', "#0072BD");
    
    % Panel Formatting
    grid on; box off;
    
    % Keep the aspect ratio square-ish and identical limits for x and y
    xlim([0, max_val * 1.1]);
    ylim([0, max_val * 1.1]);
    
    title(sprintf('%d pps', stim_rate), 'FontSize', 14);
    
    % Only add Y-label to the first panel to save space
    if r == 1
        ylabel('Spike Rate Decrement (spikes/s)', 'FontWeight', 'bold');
    end
    
    % Add X-label to the middle panel
    if r == 2
        xlabel('Onset Response Rate (spikes/s)', 'FontWeight', 'bold');
    end
    
    % Add Legend to the first panel
    if r == 1
        legend([l_ref, l_joshi, l_bruce], ...
            {"y = x", "Joshi 2017", "Bruce 2024"}, ...
            'Location', 'northwest', 'FontSize', 10);
    end
end

disp('Multi-rate evaluation complete.');