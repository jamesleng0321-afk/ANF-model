%% ========================================================================
% Heffer2010_Biophysical_Models.m
% Description: Reproduction of Heffer 2010 rate-dependent spike probability 
% changes, comparing Original Model (Joshi 2017) and Bruce (2024) Model.
% ========================================================================
clear; clc; close all;

%% 1. Load Experimental Data
load("Heffer2010.mat");                  % Ensure this file exists in your path
stim_rates = [200, 1000, 2000, 5000];    % Pulses per second (pps)
signal_length_s = 2e-3;                  % 2 ms (onset period)

%% 2. Pulse Definition
phase_len_us = 25;
ipg_us       = 8;
Fs           = 1e6; % 1 MHz sampling rate

% ---- ANF Pulse (Shared by Original and Bruce Models) ----
% Cathodic-first pulse shape mapping (25us cathodic, 8us gap, 25us anodic)
pulse_ANF = [0, -1 * ones(1, phase_len_us), ...
             zeros(1, ipg_us), ...
             1 * ones(1, phase_len_us), 0];

%% 3. Target Probability Levels
low_probs    = (0.02:0.02:0.18);
medium_probs = (0.3:0.05:0.7);
high_probs   = (0.75:0.025:0.95);

%% 4. Compute Single Pulse Thresholds 
NoiseAlpha = 0.8;
disp('Calculating Biophysical models dynamic range & thresholds...');

% Find deterministic threshold curve for a single pulse
[Level_curve, Prob_curve] = Library.FindThreshold([pulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, sigmaSingle]   = Library.FitNeuronDynamicRange(Level_curve', Prob_curve);

% Use Inverse Gaussian CDF to find exact amplitudes for target probabilities
% -------------------------------------------------------------------------
% Original ANF Levels
get_ANF_level = @(p) muSingle + sigmaSingle * sqrt(2) * erfinv(2*p - 1);
low_levels_ANF    = arrayfun(get_ANF_level, low_probs);
medium_levels_ANF = arrayfun(get_ANF_level, medium_probs);
high_levels_ANF   = arrayfun(get_ANF_level, high_probs);

% Bruce 2024 Levels (Sigma is scaled by 1/3 according to the model's RS parameter reduction)
get_Bruce_level = @(p) muSingle + (sigmaSingle / 3) * sqrt(2) * erfinv(2*p - 1);
low_levels_Bruce    = arrayfun(get_Bruce_level, low_probs);
medium_levels_Bruce = arrayfun(get_Bruce_level, medium_probs);
high_levels_Bruce   = arrayfun(get_Bruce_level, high_probs);

%% 5. Preallocate Result Structures
prob_change_ANF.low    = nan(length(stim_rates), length(low_levels_ANF));
prob_change_ANF.medium = nan(length(stim_rates), length(medium_levels_ANF));
prob_change_ANF.high   = nan(length(stim_rates), length(high_levels_ANF));

prob_change_Bruce = prob_change_ANF;

%% 6. Main Monte Carlo Simulation Loop
nTrials_MC = 50; % Monte Carlo trials (Increase to 100 for smoother curves if needed)

for s_ind = 1:length(stim_rates)
    rate = stim_rates(s_ind);
    fprintf('\nProcessing Rate: %d pps...\n', rate);
    
    % Calculate how many pulses fall into the 2ms onset window
    ipi_s = 1 / rate;
    num_pulses_in_onset = 1 + floor(signal_length_s / ipi_s);
    
    %% ---- LOW LEVELS ----
    prob_change_ANF.low(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain, pulse_ANF, low_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, low_probs);
    prob_change_Bruce.low(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain_Bruce, pulse_ANF, low_levels_Bruce, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, low_probs);
        
    %% ---- MEDIUM LEVELS ----
    prob_change_ANF.medium(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain, pulse_ANF, medium_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, medium_probs);
    prob_change_Bruce.medium(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain_Bruce, pulse_ANF, medium_levels_Bruce, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, medium_probs);
        
    %% ---- HIGH LEVELS ----
    prob_change_ANF.high(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain, pulse_ANF, high_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, high_probs);
    prob_change_Bruce.high(s_ind,:) = ...
        compute_rate_change_model(@Model_PulseTrain_Bruce, pulse_ANF, high_levels_Bruce, rate, signal_length_s, Fs, NoiseAlpha, nTrials_MC, num_pulses_in_onset, high_probs);
end

%% 7. Plotting Results
visual_off = 0.15; % Offset for visual separation of data points
fig = figure("Name","Heffer2010_Biophysical_Comparison","DefaultAxesFontSize",12, "Color","w");
fig.OuterPosition(3:4) = [1600 600];
tiledlayout(1,3, "Padding","compact");

plot_panel(data_Heffer2010(3), prob_change_ANF.high, prob_change_Bruce.high, stim_rates, visual_off, "High Probability");
plot_panel(data_Heffer2010(2), prob_change_ANF.medium, prob_change_Bruce.medium, stim_rates, visual_off, "Medium Probability");
plot_panel(data_Heffer2010(1), prob_change_ANF.low, prob_change_Bruce.low, stim_rates, visual_off, "Low Probability");

lg = legend("Orientation", "horizontal");
lg.Layout.Tile = "North";

%% ========================================================================
% Helper Functions
% ========================================================================

% Unified wrapper to compute facilitation for stochastic models via function handles
function change = compute_rate_change_model(model_func, pulse, amplitude_levels, rate, signal_len, Fs, NoiseAlpha, nTrials, num_pulses, target_probs)
    change = zeros(1, length(amplitude_levels));
    
    % Generate the normalized pulse train sequence
    Istim = Experiment.stim_PulseTrain(pulse, rate, 100, 0, signal_len, Fs);
    
    for i = 1:length(amplitude_levels)
        amp = amplitude_levels(i);
        input = Istim * amp;
        spikes_in_onset = 0;
        
        % Monte Carlo Trials
        for tr = 1:nTrials
            p_noise = Library.oneonfnoise(length(input), NoiseAlpha);
            c_noise = Library.oneonfnoise(length(input), NoiseAlpha);
            
            % Execute the passed model function
            [~, SpTimes, ~, ~] = model_func(input, p_noise, c_noise, Fs);
            
            % If at least one spike occurred within the onset window
            if any(SpTimes <= signal_len)
                spikes_in_onset = spikes_in_onset + 1;
            end
        end
        
        % Measured probability is the fraction of trials with a spike
        p_onset_measured = spikes_in_onset / nTrials;
        
        % Predicted probability based on single pulse (Assuming no facilitation)
        p_first = target_probs(i);
        predicted = min(1, p_first * num_pulses);
        
        % Facilitation = Measured - Predicted
        change(i) = p_onset_measured - predicted;
    end
end

% Function to plot data series with visual offsets
function plot_panel(data_struct, anf_change, bruce_change, rates, off, label)
    nexttile
    hold on; grid on;
    
    % 1. Experimental Data (Offset slightly to the left)
    lower_error = data_struct.change(1:4,2) - data_struct.change(5:2:end,2);
    upper_error = data_struct.change(6:2:end,2) - data_struct.change(1:4,2);
    errorbar((1:4)-off, data_struct.change(1:4,2), ...
             lower_error, upper_error, ...
             'ks','linewidth',2, 'MarkerFaceColor','k', "DisplayName", "Exp. Data (Heffer 2010)");
             
    % 2. Original ANF (Centered)
    median_a = median(anf_change,2);
    q25_a = prctile(anf_change,25,2);
    q75_a = prctile(anf_change,75,2);
    errorbar((1:4), median_a, ...
            median_a-q25_a, q75_a-median_a, ...
            'd','linewidth',2,"color","#0072BD", "DisplayName", "Original Model (ANF)");
            
    % 3. Bruce (2024) ANF (Offset slightly to the right)
    median_b = median(bruce_change,2);
    q25_b = prctile(bruce_change,25,2);
    q75_b = prctile(bruce_change,75,2);
    errorbar((1:4)+off, median_b, ...
            median_b-q25_b, q75_b-median_b, ...
            '^','linewidth',2,"color","#D95319", "DisplayName", "Bruce (2024) Model");
            
    ylabel("Spike probability change")
    xticks(1:4)
    xticklabels(string(rates))
    xlabel("Pulse rate [pps]")
    title(label)
    ylim([-0.4 1])
    xlim([0.5 4.5])
    yline(0, 'k--', 'HandleVisibility','off');
    box off; set(gca, 'TickDir', 'out');
end