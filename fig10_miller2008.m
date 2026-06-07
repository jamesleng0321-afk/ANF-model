% ========================================================================
% Script: Reproduce_and_Compare_Final.m
% Description: Reproduces Joshi et al. (2017) Figure 10 (Model Results)
% and compares it with digitized experimental data from Miller et al. (2008).
% ========================================================================

clear; clc; close all;

%% 1. Parameters Configuration
Fs = 1e6;                       % Sampling frequency (1 MHz, dt = 1 us)
NoiseAlpha = 0.8;               % 1/f noise spectral shaping parameter
stim_duration = 0.3;            % Stimulus duration (300 ms)
rates = [250, 1000, 5000];      % Stimulation rates (pulses per second)
n_trials = 50;                  % Number of Monte Carlo trials for smooth histograms
dB_above_threshold = 1;         % Stimulation level (+1 dB relative to threshold)
onset_discard_duration = 0.05;  % Discard first 50 ms to avoid onset transient effects

%% 2. Pulse Definition & Baseline Threshold Calculation
% Anodic-leading biphasic pulse: 50us anodic, 0us IPG, 50us cathodic
SinglePulse = [0, 1*ones(1,50), -1*ones(1,50), 0]; 

disp('Calculating Single Pulse Threshold (I50) using Library.FindThreshold...');
[Level, Prob] = Library.FindThreshold([SinglePulse, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, ~] = Library.FitNeuronDynamicRange(Level', Prob);
fprintf('Baseline Threshold found: %.2f uA\n', muSingle * 1e6);

%% 3. Load Miller (2008) Experimental Data from CSV
csv_file = 'Miller_2008_fig1.csv';
fprintf('Loading Experimental Data from %s...\n', csv_file);

if ~exist(csv_file, 'file')
    error('File %s not found. Please place it in the working directory.', csv_file);
end
miller_fig1 = readtable(csv_file);

%% 4. Figure Preparation
% Increased figure height from 700 to 750 to prevent title overlap
fig = figure('Position', [100, 100, 1200, 750], ...
             'Name', 'Joshi 2017 Model vs Miller 2008 Experimental Data');

%% 5. Simulation Loop for Joshi Model (Top Row)
for r_idx = 1:length(rates)
    rate = rates(r_idx);
    fprintf('Processing Model Rate: %d pps... [', rate);
    
    % Calculate exact current amplitude for +1 dB
    amplitude = muSingle * 10^(dB_above_threshold / 20);
    Istim = Experiment.stim_PulseTrain(SinglePulse, rate, 100, 0, stim_duration, Fs);
    input_current = Istim * amplitude;
    
    all_ISIs = []; 
    
    % Monte Carlo simulation loop
    for tr = 1:n_trials
        p_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
        
        [~, SpTimes, ~, ~] = Model_PulseTrain(input_current, p_noise, c_noise, Fs);
        
        SpTimes = SpTimes(:); 
        steady_SpTimes = SpTimes(SpTimes > onset_discard_duration); 
        
        % Calculate ISIs in milliseconds
        if length(steady_SpTimes) > 1
            all_ISIs = [all_ISIs; diff(steady_SpTimes) * 1000];
        end
        
        % Print progress indicator
        if mod(tr, round(n_trials/10)) == 0
            fprintf('='); 
        end
    end
    fprintf('] Done!\n');
    
    %% Plotting Model Results
    subplot(2, 3, r_idx);
    
    % Use 0.1 ms bin width to reveal sharp phase-locking peaks
    model_bin_width = 0.1;
    model_edges = 0 : model_bin_width : 20;
    model_counts = histcounts(all_ISIs, model_edges);
    
    % Normalize counts to a maximum of 1 for shape comparison
    if max(model_counts) > 0
        norm_model_counts = model_counts / max(model_counts);
    else
        norm_model_counts = model_counts;
    end
    
    histogram('BinEdges', model_edges, 'BinCounts', norm_model_counts, ...
        'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
    
    % Axis formatting
    xlim([0 20]); 
    ylim([0 1.05]);
    title(sprintf('Model: %d PPS', rate), 'FontSize', 13);
    xlabel('Inter-spike interval (ms)', 'FontSize', 11);
    
    if r_idx == 1
        ylabel('Normalized Prob.', 'FontSize', 11);
    end
    
    box off; 
    set(gca, 'TickDir', 'out', 'FontSize', 10);
    drawnow;
end

%% 6. Loop for Miller Experimental Data (Bottom Row)
fprintf('\nPlotting Digitized Experimental Data (Miller 2008)...\n');

% X-axis bins based on digitized data properties
xBins_ms = 0 : 0.05 : 30; 

for r_idx = 1:length(rates)
    subplot(2, 3, r_idx + 3); 
    rate = rates(r_idx);
    
    try
        % Extract the correct column based on the stimulation rate
        switch rate
            case 250
                data_vec = miller_fig1.x250_c2_r3;
            case 1000
                data_vec = miller_fig1.x1000_c2_r3;
            case 5000
                data_vec = miller_fig1.x5000_c2_r4;
        end
        
        % Clean up NaN values
        data_vec(isnan(data_vec)) = 0;
        
        % Normalize to a maximum of 1
        if max(data_vec) > 0
            data_vec = data_vec / max(data_vec);
        end
        
        % Plot using a filled area to mimic a continuous histogram
        fill([xBins_ms, fliplr(xBins_ms)], [data_vec', zeros(1, length(data_vec))], ...
             [0.2 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
             
        % Axis formatting
        title(sprintf('Exp (Miller): %d PPS', rate), 'FontSize', 13);
        xlabel('Inter-spike interval (ms)', 'FontSize', 11);
        
        if r_idx == 1
            ylabel('Normalized Prob.', 'FontSize', 11);
        end
        
        xlim([0 20]); 
        ylim([0 1.05]);
        
    catch ME
        % Error handling if columns are missing
        text(10, 0.5, sprintf('Data error for %d PPS', rate), ...
            'Color', 'r', 'HorizontalAlignment', 'center');
    end
    
    box off; 
    set(gca, 'TickDir', 'out', 'FontSize', 10);
end

%% 7. Add Combined Title
% Merge the main title and parameters into a two-line sgtitle to prevent overlapping
sgtitle({'Joshi 2017 Model vs Miller 2008 Experimental Data', ...
         'ISI Distributions (Normalized to Max=1) | Rate: +1 dB'}, ...
        'FontSize', 15, 'FontWeight', 'bold');