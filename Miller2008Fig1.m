%% ========================================================================
% Miller2008_Joshi_vs_Bruce_Fig1.m
% Recreate Miller 2008 Fig 1 style, comparing Joshi 2017 and Bruce 2024 models
% Uses 3 independent Figures, each structured as a 4 (currents) x 3 (time windows) grid
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% Parameter settings (Exactly matching Miller 2008 Fig 1)
% ------------------------------------------------------------------------
pulse_rates = [250, 1000, 5000];
duration_s  = 300e-3;   % 300 ms stimulus duration
nTrials     = 50;       % Monte Carlo trials (can be increased to 100+ for smoother histograms)

% Analysis time windows [Early; Mid; Late]
onset_bins = [0   12;
              4   50;
              200 300] * 1e-3; 

% Specific current levels used in Miller 2008 Fig 1 (Unit: mA)
amps_mA = {
    [1.05, 1.10, 1.15, 1.20], ... % 4 current levels for 250 pulse/s
    [1.12, 1.20, 1.30, 1.40], ... % 4 current levels for 1000 pulse/s
    [0.98, 1.05, 1.10, 1.20]  ... % 4 current levels for 5000 pulse/s
};

%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------
Fs = 1e6;
NoiseAlpha = 0.8;
% 40us positive, 0us IPG, 40us negative (Cathodic-first biphasic pulse)
SinglePulse = [0, -1*ones(1,40), 1*ones(1,40), 0]; 

% Histogram binning parameters (50 µs bin width)
bin_width = 50e-6; 
xBins_edges = 0:bin_width:22e-3; % X-axis limits 0 to 22 ms
xBins_centers = xBins_edges(1:end-1) + diff(xBins_edges)/2;

%% ========================================================================
% MAIN LOOP: Iterate over stimulus rates, current levels, and time windows
% ========================================================================
for rate_ind = 1:length(pulse_rates)
    rate = pulse_rates(rate_ind);
    amps = amps_mA{rate_ind} * 1e-3; % Convert to Amperes (A)
    
    fprintf('\nProcessing rate: %d pps...\n', rate);
    
    % Create an independent Figure for the current pulse rate
    fig = figure("Name", sprintf("Miller2008_Fig1_%d_pps", rate), "DefaultAxesFontSize", 11);
    fig.Position(3:4) = [1100 850];
    tl = tiledlayout(4, 3, "TileSpacing", "compact", "Padding", "compact");
    title(tl, sprintf('%d pulse/s (Joshi vs Bruce)', rate), 'FontSize', 16, 'FontWeight', 'bold');
    
    for level_ind = 1:length(amps)
        amp = amps(level_ind);
        fprintf('  - Applying current: %.2f mA\n', amp * 1000);
        
        Istim = Experiment.stim_PulseTrain(SinglePulse, rate, 100, 0, duration_s, Fs);
        input_current = Istim * amp;
        
        % Initialize ISI collectors for the current level
        raw_isis_Joshi = cell(1, size(onset_bins, 1));
        raw_isis_Bruce = cell(1, size(onset_bins, 1));
        
        for tr = 1:nTrials
            % Synchronize noise sequences to ensure controlled comparison
            p_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
            c_noise = Library.oneonfnoise(length(input_current), NoiseAlpha);
            
            % 1. Joshi 2017 Model
            [~, SpTimes_J, ~, ~, ~, ~] = Model_PulseTrain(input_current, p_noise, c_noise, Fs);
            if length(SpTimes_J) > 1
                SpTimes_J = SpTimes_J(:); 
                ISIs_J = diff(SpTimes_J);
                onsets_J = SpTimes_J(1:end-1);
                for b = 1:size(onset_bins, 1)
                    mask = (onsets_J >= onset_bins(b,1)) & (onsets_J <= onset_bins(b,2));
                    raw_isis_Joshi{b} = [raw_isis_Joshi{b}; ISIs_J(mask)];
                end
            end
            
            % 2. Bruce 2024 Model
            [~, SpTimes_B, ~, ~, ~, ~] = Model_PulseTrain_Bruce(input_current, p_noise, c_noise, Fs);
            if length(SpTimes_B) > 1
                SpTimes_B = SpTimes_B(:); 
                ISIs_B = diff(SpTimes_B);
                onsets_B = SpTimes_B(1:end-1);
                for b = 1:size(onset_bins, 1)
                    mask = (onsets_B >= onset_bins(b,1)) & (onsets_B <= onset_bins(b,2));
                    raw_isis_Bruce{b} = [raw_isis_Bruce{b}; ISIs_B(mask)];
                end
            end
        end

        % ---- Plot the three histograms for the current level ----
        for bin_idx = 1:size(onset_bins, 1)
            nexttile; hold on;
            
            % Calculate and normalize Joshi histogram (Orange)
            c_J = histcounts(raw_isis_Joshi{bin_idx}, xBins_edges);
            p_J = zeros(size(c_J));
            if max(c_J) > 0
                p_J = c_J / max(c_J); 
            end
            l_joshi = plot(xBins_centers*1e3, p_J, '-', 'LineWidth', 1.5, 'Color', "#D95319");
            
            % Calculate and normalize Bruce histogram (Blue)
            c_B = histcounts(raw_isis_Bruce{bin_idx}, xBins_edges);
            p_B = zeros(size(c_B));
            if max(c_B) > 0
                p_B = c_B / max(c_B); 
            end
            l_bruce = plot(xBins_centers*1e3, p_B, '-', 'LineWidth', 1.5, 'Color', "#0072BD");
            
            % Axis settings
            xlim([0 22]);
            ylim([0 1.05]);
            grid on; box off;
            
            % Annotations (Only the first row has titles)
            if level_ind == 1
                title(sprintf('%d-%d ms', onset_bins(bin_idx,1)*1000, onset_bins(bin_idx,2)*1000));
            end
            
            % Annotations (Only the last row has X-axis labels)
            if level_ind == 4
                xlabel("Inter-spike interval (ms)");
            end
            
            % Annotations (Only the first column has Y-axis labels and current levels)
            if bin_idx == 1
                ylabel(sprintf('%.2f mA\nNorm. Count', amp*1000), 'FontWeight', 'bold');
            end
        end
    end
    
    % Add legend at the top of each Figure
    lg = legend([l_joshi, l_bruce], ["Joshi 2017", "Bruce 2024"], "Orientation", "horizontal");
    lg.Layout.Tile = "north";
end

disp('Plotting complete.');