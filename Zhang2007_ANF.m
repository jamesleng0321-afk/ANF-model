%% ========================================================================
% Model Comparison aPSTH Script (Strictly Zhang 2007 Conditions)
% Description: Evaluates and compares the aPSTH responses of the Original 
%              Model (Joshi 2017) and the Bruce (2024) Improved Model.
%% ========================================================================
clear; clc; close all;

%% 1. Parameter Settings
Fs = 1e6;                       % Sampling frequency (1 MHz)
NoiseAlpha = 0.8;               % 1/f noise spectral shaping parameter
stim_duration = 0.3;            % Stimulus duration (seconds)
nTrials = 100;                  % Number of Monte Carlo trials (50 or 100 recommended for smooth curves)

% Time bin definitions for aPSTH calculation (Based on Zhang 2007)
x_bins = [0;4;12;24;36;48;100;200;300] * 1e-3; 
x_values = x_bins(2:end) - diff(x_bins)/2; % Bin centers for plotting

% Experimental conditions: Stimulation rates (pps) and levels (dB re I50)
pulse_rates = [250, 1000, 5000];
amplitudes_dB = [0.5, 0.5, 3.07;   
                 1.3, 6.8, 5.96];  

% Preallocate matrices for storing PSTH results (Rates, Amplitudes, Bins)
mean_rate_orig = nan(length(pulse_rates), 2, length(x_bins)-1);
mean_rate_bruce = nan(length(pulse_rates), 2, length(x_bins)-1);

%% 2. Accurately Calculate Baseline Single Pulse Threshold (I50)
% [Core Fix]: Restore the 40us cathodic-first pulse with an 8us interphase gap from Zhang2007.m
SinglePulse = [0, -1*ones(1,40e-6*Fs), zeros(1,8e-6*Fs), 1*ones(1,40e-6*Fs), 0];

disp('Calculating Single Pulse Threshold (I50) using Library.FindThreshold...');
[Level, Prob] = Library.FindThreshold([SinglePulse, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, ~] = Library.FitNeuronDynamicRange(Level', Prob);
I50 = muSingle;
fprintf('Baseline Threshold (I50) found: %.4f uA\n', I50 * 1e6);

%% 3. Monte Carlo Main Loop
for p_ind = 1:length(pulse_rates)
    rate = pulse_rates(p_ind);
    
    for amp_ind = 1:2
        % Calculate absolute current amplitude for the current dB level
        amp_dB_curr = amplitudes_dB(amp_ind, p_ind);
        amplitude = I50 * 10^(amp_dB_curr / 20);
        
        fprintf('Processing Rate: %d pps, Amp: %.2f dB re I50...\n', rate, amp_dB_curr);
        
        % Construct the continuous pulse train stimulus
        Istim_base = Experiment.stim_PulseTrain(SinglePulse, rate, 100, 0, stim_duration, Fs);
        input_current = Istim_base * amplitude; 
        
        allSpTimes_orig = [];
        allSpTimes_bruce = [];
        
        for i = 1:nTrials
            % Share identical noise sequences within the same Trial for fair control
            p_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
            c_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
            
            % Run the Original Model (corresponds to ANF in the reference plot)
            [~, SpTimes_orig_n, ~, ~] = Model_PulseTrain(input_current, p_noise, c_noise, Fs);
            allSpTimes_orig = [allSpTimes_orig; SpTimes_orig_n(:)]; 
            
            % Run the Bruce (2024) Improved Model
            [~, SpTimes_bruce_n, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_current, p_noise, c_noise, Fs);
            allSpTimes_bruce = [allSpTimes_bruce; SpTimes_bruce_n(:)];
        end
        
        % Calculate binned firing rates with strict dimension control to prevent flat-line bugs
        [mean_rate_orig(p_ind, amp_ind, :)] = compute_binned_rate(allSpTimes_orig, x_bins, nTrials);
        [mean_rate_bruce(p_ind, amp_ind, :)] = compute_binned_rate(allSpTimes_bruce, x_bins, nTrials);
    end
end

%% 4. Plot Comparison Chart
% Attempt to load digitized experimental data from CSV
try
    data = readtable("Zhang_2007_fig2.csv");
    data_inds = [4, 7, 9; 2, 5, 8]; 
    has_historical = true;
catch
    warning('Zhang_2007_fig2.csv not found. Historical data line will be omitted.');
    has_historical = false;
end

fig = figure("Name","Neuron_Model_Comparison_MonteCarlo", "Color","w");
fig.OuterPosition(3:4) = [1200, 700];
t = tiledlayout(2,3,"Padding","none");
t.TileIndexing = 'columnmajor';

for p_ind = 1:length(pulse_rates)
    for amp_ind = 1:2
        nexttile;
        hold on;
        
        y_orig = squeeze(mean_rate_orig(p_ind, amp_ind, :));
        y_bruce = squeeze(mean_rate_bruce(p_ind, amp_ind, :));
        
        % Plotting curves: Original (ANF) with blue circles, Bruce with red squares
        plot(x_values, y_orig, 'o-', 'LineWidth',2, "Color","#0072BD", "DisplayName", "Original Model (ANF)");
        plot(x_values, y_bruce, 's-', 'LineWidth',2, "Color","#D95319", "DisplayName", "Bruce (2024) Model");
        
        % Enforce inclusion of the historical data markers
        if has_historical
            y_data = data{:, data_inds(amp_ind, p_ind)};
            plot(x_values, y_data, 'xk', 'LineWidth',2, 'MarkerSize',8, "DisplayName", "Historical Data");
            max_y = max([y_orig; y_bruce; y_data; 100]);
        else
            max_y = max([y_orig; y_bruce; 100]);
        end
        
        xlabel("Time [s]");
        ylabel("Spike rate [spikes/s]");
        title(sprintf("%d pps, %.1f dB re I_{50}", pulse_rates(p_ind), amplitudes_dB(amp_ind, p_ind)));
        grid on;
        ylim([0, max_y * 1.1]);
    end
end

% Create a unified legend at the top of the tiled layout
leg = legend(gca, "Orientation","horizontal");
leg.Layout.Tile = 'north';

%% ========================================================================
% Helper Function: Binned Rate Calculator with Strict Vector Dimensions
% ========================================================================
function binned_rate = compute_binned_rate(all_sp_times, x_bins, n_trials)
    [counts, ~] = histcounts(all_sp_times, x_bins);
    
    % Force conversion to column vectors to prevent matrix broadcasting behavior
    counts = counts(:); 
    x_bins = x_bins(:);
    win_width = diff(x_bins);
    
    % Strictly execute element-wise division for accurate rate calculation
    binned_rate = counts ./ (n_trials .* win_width); 
end