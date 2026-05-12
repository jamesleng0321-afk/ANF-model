%% ========================================================================
% Main Comparison Script: psBLIF vs ANF
% Normalized comparison of probabilistic and stochastic auditory models
% ========================================================================
clear; clc; close all;

%% 1. Parameters Setup
% Data and Bin definitions
data = readtable("Zhang_2007_fig2.csv");
data_inds = [4, 7, 9; 2, 5, 8]; 
x_bins = [0;4;12;24;36;48;100;200;300] * 1e-3; % seconds
x_values = x_bins(2:end) - diff(x_bins)/2;
pulse_rates = [250, 1000, 5000];

% Amplitudes in dB relative to threshold

amplitudes_dB = [0.5, 1.3, 3.07; 
                 0.5, 1.3, 5.96]; 

% Preallocate matrices for Spike Rates

mean_spiking_probability_psBLIF = nan(length(pulse_rates), 2, length(x_bins)-1);
mean_spiking_probability_ANF    = nan(length(pulse_rates), 2, length(x_bins)-1);

%% 2. Thresholding
% ---- psBLIF Threshold ----
C = psBLIF_default_parameters();
pulse_psBLIF.pulse_onset = 0; pulse_psBLIF.positive_duration = 40e-6; 
pulse_psBLIF.positive_amplitude = 1; pulse_psBLIF.interphase_gap = 0; 
pulse_psBLIF.negative_duration = 40e-6; pulse_psBLIF.negative_amplitude = -1;

eval_prob = @(amp) psBLIF_final_pulse_prob( ...
                    @(a) psBLIF_wrapper_scale_all(a, pulse_psBLIF, C), amp);
psBLIF_I50 = psBLIF_find_threshold(eval_prob, 0.5);

% ---- ANF Threshold ----
% Using the logic from your provided script

Fs = 1e6;
NoiseAlpha = 0.8;
SinglePulse = [0, -1*ones(1,40e-6*Fs), zeros(1,8e-6*Fs), 1*ones(1,40e-6*Fs), 0];
[Level,Probability] = Library.FindThreshold([SinglePulse, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, sigmaSingle] = Library.FitNeuronDynamicRange(Level', Probability);
ANF_I50 = muSingle; % [EN] ANF baseline threshold / [ZH] ANF 基准阈值

%% 3. MAIN LOOP
% Loop through rates and amplitudes to get responses for both models

nTrials_ANF = 10; % Number of Monte Carlo trials 

for p_ind = 1:length(pulse_rates)
    rate = pulse_rates(p_ind);
    fprintf('\nProcessing Rate: %d pps...\n', rate); 
    
    % Build psBLIF train

    pulse_train = build_psBLIF_train(0.3, rate, pulse_psBLIF);
    
    for amp_ind = 1:2
        % ---- Calculate Absolute Amplitudes
        amp_dB_curr = amplitudes_dB(amp_ind, p_ind);
        amp_psBLIF = psBLIF_I50 * 10^(amp_dB_curr / 20);
        amp_ANF    = ANF_I50 * 10^(amp_dB_curr / 20);
        
        % ========================================================
        % Model 1: psBLIF (Analytical)
        % ========================================================
        [~, out_ps] = psBLIF_wrapper_scale_all(amp_psBLIF, pulse_train, C);
        t_max = pulse_train(end).pulse_onset + 1e-2;
        dist_psBLIF = psBLIF_get_spike_distribution_fast(out_ps, pulse_train, t_max, 1e6);
        
        for ind = 1:(length(x_bins)-1)
            i1 = round(x_bins(ind) * 1e6) + 1;
            i2 = round(x_bins(ind + 1) * 1e6);
            mean_spiking_probability_psBLIF(p_ind, amp_ind, ind) = ...
                sum(dist_psBLIF(i1:i2)) / (x_bins(ind+1) - x_bins(ind));
        end
        
        % ========================================================
        % Model 2: ANF (Stochastic)
        % ========================================================
        % Call our custom wrapper function to get binned aPSTH rates

        stim_duration = 0.3;
        rate_ANF_binned = get_ANF_aPSTH_rates_Wrapper(rate, amp_ANF, SinglePulse, ...
                                    stim_duration, Fs, NoiseAlpha, nTrials_ANF, x_bins);
        
        mean_spiking_probability_ANF(p_ind, amp_ind, :) = rate_ANF_binned;
    end
end

%% 4. Plotting (Normalized Output) 
fig = figure("Name","Model_Comparison_Normalized", "DefaultAxesFontSize",13);
fig.OuterPosition(3:4) = [1200, 700];
t = tiledlayout(2,3,"Padding","none");
t.TileIndexing = 'columnmajor';

for p_ind = 1:length(pulse_rates)
    for amp_ind = 1:2
        ax = nexttile;
        hold on;
        
        % Extract data for current condition
        y_psBLIF = squeeze(mean_spiking_probability_psBLIF(p_ind, amp_ind, :));
        y_ANF    = squeeze(mean_spiking_probability_ANF(p_ind, amp_ind, :));
        y_data   = data{:, data_inds(amp_ind, p_ind)};
        
        % ========================================================
        % OUTPUT NORMALIZATION (Normalize to Max Rate)
        % ========================================================
        y_psBLIF_norm = y_psBLIF / max(y_psBLIF);
        y_ANF_norm    = y_ANF / max(y_ANF);
        y_data_norm   = y_data / max(y_data);
        
        % Plotting the normalized curves
        plot(x_values, y_psBLIF_norm, 'o-', 'LineWidth',2, "Color","#D95319", "DisplayName", "psBLIF");
        plot(x_values, y_ANF_norm, 's-', 'LineWidth',2, "Color","#0072BD", "DisplayName", "ANF");
        plot(x_values, y_data_norm, 'xk', 'LineWidth',2, 'MarkerSize',8, "DisplayName", "Data");
        
        xlabel("Time [s]");
        ylabel("Normalized Spike rate");
        title(sprintf("%d pps, %.1f dB re I_{50}", pulse_rates(p_ind), amplitudes_dB(amp_ind, p_ind)));
        grid on;
        ylim([0, 1.1]); % [EN] Y-axis 0 to 1.1 since it's normalized / [ZH] 因为归一化了，Y轴范围为 0 到 1.1
    end
end
leg = legend(ax, "Orientation","horizontal");
leg.Layout.Tile = 'north';


%% ========================================================================
% Helper Functions
% ========================================================================

% Wrapper Function for ANF Model
function aPSTH_rates = get_ANF_aPSTH_rates_Wrapper(stim_rate, amplitude, SinglePulse, stim_duration, Fs, NoiseAlpha, nTrials, x_bins)
    % Generate the pulse train with the absolute amplitude
    Istim = Experiment.stim_PulseTrain(SinglePulse, stim_rate, 100, 0, stim_duration, Fs);
    input = Istim * amplitude; 
    
    allSpikeTimes = [];
    
    % Run Monte Carlo trials

    for i = 1:nTrials
        p_noise = Library.oneonfnoise(length(input), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(input), NoiseAlpha);
        
        % Call your local model function
        [~, SpTimes, ~, ~] = Model_PulseTrain(input, p_noise, c_noise, Fs);
        allSpikeTimes = [allSpikeTimes; SpTimes(:)]; 
    end
    
    % Calculate aPSTH based on the exact same x_bins as psBLIF
    aPSTH_rates = zeros(1, length(x_bins)-1);
    
    for i = 1:(length(x_bins)-1)
        t_start = x_bins(i);
        t_end   = x_bins(i+1);
        % Find spikes within current bin
        spike_count_in_win = sum(allSpikeTimes >= t_start & allSpikeTimes < t_end);
        win_width = t_end - t_start;
        % Convert count to rate (Spikes/second)
        aPSTH_rates(i) = spike_count_in_win / (nTrials * win_width);
    end
end

% ... (Keep your other helper functions like build_psBLIF_train, psBLIF_wrapper_scale_all here)


%% ========================================================================
% Helper Functions
% ========================================================================

function p = psBLIF_final_pulse_prob(fun, amp)
[~, out] = fun(amp);
p = sum([out{end}.path_prob]);
end

function pulse_train = build_psBLIF_train(duration_s, rate, pulse)
ipi = 1/rate;
n   = floor(duration_s / ipi);
pulse_train = repmat(pulse, 1, n);
for k = 1:n
    pulse_train(k).pulse_onset = (k-1)*ipi;
end
end

function [history, out] = psBLIF_wrapper_scale_all(amplitude, pulse_train, C)
scaled = pulse_train;
for k = 1:numel(scaled)
    scaled(k).positive_amplitude = scaled(k).positive_amplitude * amplitude;
    scaled(k).negative_amplitude = scaled(k).negative_amplitude * amplitude;
end
[history, out] = psBLIF(scaled, C);
end

function rate = psBLIF_spike_rate(amp, pulse_train, C, bin_width)
[history, alifp_ret] = psBLIF_wrapper_scale_all(amp, pulse_train, C);
p_vec = psBLIF_per_pulse_prob(alifp_ret);
rate = sum(p_vec)/bin_width;
end