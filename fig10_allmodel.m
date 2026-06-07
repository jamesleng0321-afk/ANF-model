% ========================================================================
% Script: Compare_All_Models.m
% Description: Comprehensive comparison of ISI Histograms.
% Row 1: Original Joshi 2017 Model
% Row 2: Improved Bruce 2024 Model (Scaled RS + Random Init)
% Row 3: Miller 2008 Experimental Data
% ========================================================================

clear; clc; close all;

%% 1. Parameters Configuration
Fs = 1e6;                       % Sampling frequency (1 MHz, dt = 1 us)
NoiseAlpha = 0.8;               % 1/f noise spectral shaping parameter
stim_duration = 0.3;            % Stimulus duration (300 ms)
rates = [250, 1000, 5000];      % Stimulation rates (pulses per second)
n_trials = 50;                  % Number of Monte Carlo trials
dB_above_threshold = 1;         % Stimulation level (+1 dB relative to threshold)
onset_discard_duration = 0.05;  % Discard first 50 ms (onset transient)

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
% Create a larger figure (3x3 grid)
fig = figure('Position', [50, 50, 1200, 900], ...
             'Name', 'Comprehensive Comparison: Joshi 2017 vs Bruce 2024 vs Miller 2008');

%% 5. Main Simulation Loop for Both Models
for r_idx = 1:length(rates)
    rate = rates(r_idx);
    fprintf('\n========================================\n');
    fprintf('Processing Rate: %d pps...\n', rate);
    
    % Calculate exact current amplitude for +1 dB
    amplitude = muSingle * 10^(dB_above_threshold / 20);
    Istim = Experiment.stim_PulseTrain(SinglePulse, rate, 100, 0, stim_duration, Fs);
    input_current = Istim * amplitude;
    
    all_ISIs_orig = []; 
    all_ISIs_impr = [];
    
    fprintf('Simulating Original & Improved Models [');
    for tr = 1:n_trials
        p_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
        
        % --- Run Original Model (Joshi 2017) ---
        [~, SpTimes_orig, ~, ~] = Model_PulseTrain(input_current, p_noise, c_noise, Fs);
        SpTimes_orig = SpTimes_orig(:);
        steady_orig = SpTimes_orig(SpTimes_orig > onset_discard_duration); 
        if length(steady_orig) > 1
            all_ISIs_orig = [all_ISIs_orig; diff(steady_orig) * 1000];
        end
        
        % --- Run Improved Model (Bruce 2024) ---
        [~, SpTimes_impr, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_current, p_noise, c_noise, Fs);
        SpTimes_impr = SpTimes_impr(:);
        steady_impr = SpTimes_impr(SpTimes_impr > onset_discard_duration);
        if length(steady_impr) > 1
            all_ISIs_impr = [all_ISIs_impr; diff(steady_impr) * 1000];
        end
        
        if mod(tr, round(n_trials/10)) == 0, fprintf('='); end
    end
    fprintf('] Done!\n');
    
    % Common histogram edges (0.1 ms resolution)
    edges = 0 : 0.1 : 20;
    
    %% Plot Row 1: Original Model
    subplot(3, 3, r_idx);
    counts_orig = histcounts(all_ISIs_orig, edges);
    norm_orig = counts_orig / max([1, max(counts_orig)]); % Normalize to max=1
    
    histogram('BinEdges', edges, 'BinCounts', norm_orig, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    xlim([0 20]); ylim([0 1.05]);
    title(sprintf('Original Model: %d PPS', rate), 'FontSize', 12);
    if r_idx == 1, ylabel('Norm. Prob.', 'FontSize', 11, 'FontWeight', 'bold'); end
    box off; set(gca, 'TickDir', 'out');
    
    %% Plot Row 2: Improved Model (Bruce 2024)
    subplot(3, 3, r_idx + 3);
    counts_impr = histcounts(all_ISIs_impr, edges);
    norm_impr = counts_impr / max([1, max(counts_impr)]); % Normalize to max=1
    
    histogram('BinEdges', edges, 'BinCounts', norm_impr, 'FaceColor', [0.8 0.5 0.4], 'EdgeColor', 'none');
    xlim([0 20]); ylim([0 1.05]);
    title(sprintf('Improved Model: %d PPS', rate), 'FontSize', 12);
    if r_idx == 1, ylabel('Norm. Prob.', 'FontSize', 11, 'FontWeight', 'bold'); end
    box off; set(gca, 'TickDir', 'out');
    
    drawnow;
end

%% 6. Plot Row 3: Miller 2008 Experimental Data
fprintf('\nPlotting Digitized Experimental Data (Miller 2008)...\n');
xBins_ms = 0 : 0.05 : 30; 

for r_idx = 1:length(rates)
    subplot(3, 3, r_idx + 6); 
    rate = rates(r_idx);
    
    try
        switch rate
            case 250
                data_vec = miller_fig1.x250_c2_r3;
            case 1000
                data_vec = miller_fig1.x1000_c2_r3;
            case 5000
                data_vec = miller_fig1.x5000_c2_r4;
        end
        data_vec(isnan(data_vec)) = 0;
        if max(data_vec) > 0, data_vec = data_vec / max(data_vec); end
        
        fill([xBins_ms, fliplr(xBins_ms)], [data_vec', zeros(1, length(data_vec))], ...
             [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
             
        title(sprintf('Exp (Miller): %d PPS', rate), 'FontSize', 12);
        xlabel('Inter-spike interval (ms)', 'FontSize', 11);
        if r_idx == 1, ylabel('Norm. Prob.', 'FontSize', 11, 'FontWeight', 'bold'); end
        xlim([0 20]); ylim([0 1.05]);
        
    catch
        text(10, 0.5, 'Data Error', 'Color', 'r', 'HorizontalAlignment', 'center');
    end
    box off; set(gca, 'TickDir', 'out');
end

%% 7. Add Combined Title
sgtitle({'ISI Distributions: Original vs. Improved Model vs. Experimental Data', ...
         '(Normalized Probabilities | Rate: +1 dB)'}, ...
        'FontSize', 16, 'FontWeight', 'bold');