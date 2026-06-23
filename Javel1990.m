%% ========================================================================
% Evaluate_Javel1990_BouletStyle.m
% Reproduction of Javel 1990: Mean firing probability vs stimulus level
% Comparison: Joshi 2017 vs Bruce 2024
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Load Experimental Data (Javel 1990)
% ------------------------------------------------------------------------
try
    load('Javel1990.mat');
    has_javel_data = true;
    fprintf('Loaded Javel1990.mat successfully.\n');
catch
    warning('Javel1990.mat not found. Will plot models only.');
    has_javel_data = false;
end

colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];

%% ------------------------------------------------------------------------
% 2. Stimulation Parameters
% ------------------------------------------------------------------------
pulse_rates = [100, 200, 400, 800];   % pps (Pulses per second)
duration_s  = 0.1;                    % 100 ms stimulus duration
nTrials     = 15;                     % Monte Carlo trials

Fs = 1e6;
NoiseAlpha = 0.8;

% Javel 1990 used 50 µs/phase biphasic pulses
% Define: 50us negative, 0us IPG, 50us positive
SinglePulse = [0, -1*ones(1,50), 1*ones(1,50), 0];

% Define current sweep range: 0.4 mA to 2.0 mA
amplitudes_mA = 0.4 : 0.05 : 2.0;
amplitudes_A  = amplitudes_mA * 1e-3; % Convert to Amperes

%% ------------------------------------------------------------------------
% 3. Preallocate Storage
% ------------------------------------------------------------------------
mean_prob_Joshi = zeros(length(pulse_rates), length(amplitudes_A));
mean_prob_Bruce = zeros(length(pulse_rates), length(amplitudes_A));

%% ========================================================================
% MAIN LOOP: Sweep over frequencies and amplitudes
% ========================================================================
fprintf('=== Simulating Mean Firing Probabilities ===\n');

for p_ind = 1:length(pulse_rates)
    rate = pulse_rates(p_ind);
    fprintf('\nProcessing Pulse Rate: %d pps...\n', rate);
    
    % Number of total pulses in one stimulus presentation
    num_pulses = floor(duration_s * rate);
    
    % Build normalized stimulus train
    Istim_base = Experiment.stim_PulseTrain(SinglePulse, rate, 100, 0, duration_s, Fs);

    for a_ind = 1:length(amplitudes_A)
        input_I = Istim_base * amplitudes_A(a_ind);
        
        spikes_J = 0;
        spikes_B = 0;
        
        for tr = 1:nTrials
            % Ensure both models receive the same random noise trace per trial
            p_n = Library.oneonfnoise(length(input_I), NoiseAlpha);
            c_n = Library.oneonfnoise(length(input_I), NoiseAlpha);
            
            % Run Joshi 2017
            [~, SpT_J, ~, ~, ~, ~] = Model_PulseTrain(input_I, p_n, c_n, Fs);
            spikes_J = spikes_J + length(SpT_J);
            
            % Run Bruce 2024
            [~, SpT_B, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_I, p_n, c_n, Fs);
            spikes_B = spikes_B + length(SpT_B);
        end
        
        % Mean probability per pulse = Total Spikes / (Total Pulses * Trials)
        % Note: Probability can exceed 1 if multi-spiking occurs, cap it at 1 for probability curve
        mean_prob_Joshi(p_ind, a_ind) = min(1.0, spikes_J / (num_pulses * nTrials));
        mean_prob_Bruce(p_ind, a_ind) = min(1.0, spikes_B / (num_pulses * nTrials));
    end
end

%% ------------------------------------------------------------------------
% 4. Calculate Reference Thresholds (I50 at 100 pps)
% ------------------------------------------------------------------------
% We use the 100 pps condition to find the 50% probability threshold
I50_Joshi = find_50_threshold(amplitudes_A, mean_prob_Joshi(1,:));
I50_Bruce = find_50_threshold(amplitudes_A, mean_prob_Bruce(1,:));

fprintf('\nReference Threshold (100 pps):\n');
fprintf('  Joshi 2017: %.3f mA\n', I50_Joshi * 1000);
fprintf('  Bruce 2024: %.3f mA\n', I50_Bruce * 1000);

%% ========================================================================
% 5. Plotting (Recreating Javel 1990 Style)
% ========================================================================
fig = figure("Name", "Javel1990_Joshi_vs_Bruce", "DefaultAxesFontSize", 13);
fig.InnerPosition(3:4) = [800 600];
hold on;

lgh = []; % Array to store legend handles
lg_labels = {}; % Array to store legend labels

for p_ind = 1:length(pulse_rates)
    rate = pulse_rates(p_ind);
    
    % ---- A. Plot Experimental Data ----
    if has_javel_data
        [data_amp_db, data_prob] = extract_javel_data(data_Javel1990, rate);
        data_amp_linear = 10.^(data_amp_db/20);
        
        % Calculate data threshold relative to 100 pps
        if rate == 100
            try
                % Assuming aLIFP_fit_gaussian exists in user's original workspace
                [data_threshold, ~] = aLIFP_fit_gaussian(data_amp_linear, data_prob);
            catch
                % Fallback: manual 50% cross
                data_threshold = find_50_threshold(data_amp_linear', data_prob');
            end
        end
        
        l_exp = plot(20*log10(data_amp_linear / data_threshold), data_prob, ...
             'x', 'Color', colors(p_ind), 'MarkerSize', 8, 'LineWidth', 2);
    end

    % ---- B. Plot Joshi 2017 ----
    l_j = plot(20*log10(amplitudes_A / I50_Joshi), mean_prob_Joshi(p_ind,:), ...
         '--', 'Color', colors(p_ind), 'LineWidth', 2);

    % ---- C. Plot Bruce 2024 ----
    l_b = plot(20*log10(amplitudes_A / I50_Bruce), mean_prob_Bruce(p_ind,:), ...
         '-', 'Color', colors(p_ind), 'LineWidth', 2);
         
    % Collect handles for legend
    if has_javel_data
        lgh = [lgh, l_exp, l_j, l_b];
        lg_labels = [lg_labels, ...
            sprintf("Data %d pps", rate), ...
            sprintf("Joshi %d pps", rate), ...
            sprintf("Bruce %d pps", rate)];
    else
        lgh = [lgh, l_j, l_b];
        lg_labels = [lg_labels, ...
            sprintf("Joshi %d pps", rate), ...
            sprintf("Bruce %d pps", rate)];
    end
end

% Formatting
ylabel("Mean Firing Probability", 'FontWeight', 'bold');
xlabel("Stimulus Level [dB re 100 pps threshold]", 'FontWeight', 'bold');
title('Mean Firing Probability vs Stimulus Level (Javel 1990)');
xlim([-2 8]);
ylim([0 1.05]);
grid on; box off;

legend(lgh, lg_labels, "Location", "northoutside", "Orientation", "horizontal", "NumColumns", 3, 'FontSize', 10);

disp('Plotting complete.');

%% ========================================================================
% Helper Functions
% ========================================================================

function threshold = find_50_threshold(amps, probs)
    % A robust helper to find the exact amplitude that crosses 0.5 probability
    [unique_probs, idx] = unique(probs);
    unique_amps = amps(idx);
    
    if max(unique_probs) < 0.5
        warning('Model did not reach 50% probability in the given amplitude range.');
        threshold = NaN;
    else
        threshold = interp1(unique_probs, unique_amps, 0.5, 'linear', 'extrap');
    end
end

function [amp, prob] = extract_javel_data(data_struct, rate)
    % Extractor from original script for the Javel dataset structure
    switch rate
        case 100
            amp  = data_struct.hundred(:,1);
            prob = data_struct.hundred(:,2);
        case 200
            amp  = data_struct.twoH(:,1);
            prob = data_struct.twoH(:,2);
        case 400
            amp  = data_struct.fourH(:,1);
            prob = data_struct.fourH(:,2);
        case 800
            amp  = data_struct.eightH(:,1);
            prob = data_struct.eightH(:,2);
        otherwise
            error("Unsupported pulse rate.")
    end
end