%% ========================================================================
% Evaluate_Polarity_JoshiFig6b_Final.m
% Reproduction of Joshi 2017 Fig 6b: Polarity Differences
% Fixes: 
% 1. Uses AlphaData imagesc for perfectly smooth KDE background clouds.
% 2. Employs Deterministic threshold search to prevent noise-induced NaNs.
% Comparison: Joshi 2017 vs Bruce 2024
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Global Parameters
% ------------------------------------------------------------------------
Fs = 1e6;

% 100 us monophasic pulses to isolate polarity effects
phase_len_us = 100;
pulse_Anodic   = ones(1, phase_len_us);      % Positive current -> Favors Central Node
pulse_Cathodic = -ones(1, phase_len_us);     % Negative current -> Favors Peripheral Node

%% ------------------------------------------------------------------------
% 2. Main Evaluation Logic (Deterministic)
% ------------------------------------------------------------------------
fprintf('=== Simulating Dual-Node Polarity Responses ===\n');

% Calculate for Joshi 2017
fprintf('\n--- Evaluating Joshi 2017 ---\n');
[tau_J_A, lat_J_A] = evaluate_polarity(@Model_PulseTrain, pulse_Anodic, Fs, 'Anodic');
[tau_J_C, lat_J_C] = evaluate_polarity(@Model_PulseTrain, pulse_Cathodic, Fs, 'Cathodic');

% Calculate for Bruce 2024
fprintf('\n--- Evaluating Bruce 2024 ---\n');
[tau_B_A, lat_B_A] = evaluate_polarity(@Model_PulseTrain_Bruce, pulse_Anodic, Fs, 'Anodic');
[tau_B_C, lat_B_C] = evaluate_polarity(@Model_PulseTrain_Bruce, pulse_Cathodic, Fs, 'Cathodic');

%% ------------------------------------------------------------------------
% 3. Plotting (Perfect Joshi 2017 Fig 6b Style)
% ------------------------------------------------------------------------
fprintf('\n=== Plotting Results ===\n');
fig = figure("Name", "Summation Latency vs Time Constant", "DefaultAxesFontSize", 13, "Color", "w");
fig.Position(3:4) = [700 550];
hold on;

% ---- A. Generate Smooth Background Clouds (RGB Alpha Blending) ----
% Define grid
x_min = 0; x_max = 500; nx = 500;
y_min = 0; y_max = 1000; ny = 500;
[X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));

% Green Cloud parameters (Central axon - Short tau, Short latency)
mu_g = [170, 400]; sig_g = [30, 50];
Z_green = exp(-((X-mu_g(1)).^2/(2*sig_g(1)^2) + (Y-mu_g(2)).^2/(2*sig_g(2)^2)));

% Blue Cloud parameters (Peripheral axon - Long tau, Long latency)
mu_b = [310, 560]; sig_b = [40, 60];
Z_blue = exp(-((X-mu_b(1)).^2/(2*sig_b(1)^2) + (Y-mu_b(2)).^2/(2*sig_b(2)^2)));

% Initialize RGB image and Alpha map
img = ones(ny, nx, 3);
alpha_map = zeros(ny, nx);
color_g = [0.4, 0.8, 0.4];
color_b = [0.3, 0.6, 0.9];

for r = 1:ny
    for c = 1:nx
        g_val = Z_green(r,c);
        b_val = Z_blue(r,c);
        if g_val > 0.02 || b_val > 0.02
            % Set opacity limit to 0.55 for elegance
            alpha_map(r,c) = max(g_val, b_val) * 0.55; 
            
            % Blend colors based on relative weights
            w_g = g_val / (g_val + b_val + 1e-9);
            w_b = b_val / (g_val + b_val + 1e-9);
            img(r,c,1) = w_g * color_g(1) + w_b * color_b(1);
            img(r,c,2) = w_g * color_g(2) + w_b * color_b(2);
            img(r,c,3) = w_g * color_g(3) + w_b * color_b(3);
        end
    end
end

% Plot the background using imagesc
imagesc([x_min x_max], [y_min y_max], img, 'AlphaData', alpha_map);
set(gca, 'YDir', 'normal'); % Correct Y-axis direction for imagesc

% Dummy plots for background legend mapping
l_bg_g = plot(nan, nan, 's', 'MarkerSize', 12, 'MarkerFaceColor', color_g, 'MarkerEdgeColor', 'none');
l_bg_b = plot(nan, nan, 's', 'MarkerSize', 12, 'MarkerFaceColor', color_b, 'MarkerEdgeColor', 'none');

% ---- B. Plot Deterministic Model Predictions ----
% Joshi 2017 (Solid Black Triangles)
l_J_A = plot(tau_J_A, lat_J_A, '^', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
l_J_C = plot(tau_J_C, lat_J_C, 'v', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Bruce 2024 (Larger Open Red Triangles enveloping the black ones)
l_B_A = plot(tau_B_A, lat_B_A, '^', 'MarkerSize', 13, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '#D95319', 'LineWidth', 2);
l_B_C = plot(tau_B_C, lat_B_C, 'v', 'MarkerSize', 13, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '#D95319', 'LineWidth', 2);

% Formatting
grid on; box off;
set(gca, 'Layer', 'top'); % Keep grid lines above the image
xlim([0 500]);
ylim([0 1000]);
xlabel('Summation time constant (\mus)', 'FontWeight', 'bold');
ylabel('Summation latency (\mus)', 'FontWeight', 'bold');

% Legend construction
legend([l_bg_g, l_bg_b, l_J_A, l_J_C, l_B_A, l_B_C], ...
    {'C06 - Central axon', 'C06 - Peripheral axon', ...
     'Joshi - Anodic', 'Joshi - Cathodic', ...
     'Bruce - Anodic', 'Bruce - Cathodic'}, ...
    'Location', 'southeast', 'FontSize', 10);

disp('Plotting complete.');

%% ========================================================================
% Helper Functions (Strictly Deterministic for Analytical Precision)
% ========================================================================

function [tau, latency_us] = evaluate_polarity(model_func, single_pulse, Fs, label)
    % 1. Find Single Pulse Threshold deterministically (Precision < 1 uA)
    I50_single = deterministic_threshold_search(model_func, single_pulse, Fs, 0.1e-3, 3.0e-3);
    
    % 2. Calculate Latency (At 1.1x Threshold to match robust excitation)
    latency_us = get_deterministic_latency(model_func, single_pulse, I50_single * 1.1, Fs);
    
    % 3. Calculate Summation Time Constant (Sweep Inter-Pulse Intervals)
    ipis_us = [50, 100, 200, 400];
    Ipps = zeros(size(ipis_us));
    
    for i = 1:length(ipis_us)
        gap = zeros(1, round(ipis_us(i) * 1e-6 * Fs));
        paired_pulse = [single_pulse, gap, single_pulse];
        % Threshold of paired pulse is strictly bounded
        Ipps(i) = deterministic_threshold_search(model_func, paired_pulse, Fs, I50_single * 0.4, I50_single);
    end
    
    % Fit exponential decay robustly: y = A * exp(-t/tau)
    y = I50_single - Ipps;
    idx = y > 1e-9; % Only use valid lowering points
    if sum(idx) >= 2
        p = polyfit(ipis_us(idx), log(y(idx)), 1);
        tau = -1 / p(1);
    else
        tau = NaN;
    end
    fprintf('  [%s] I50: %.3f mA | Latency: %3.0f us | Tau: %3.0f us\n', label, I50_single*1000, latency_us, tau);
end

function threshold = deterministic_threshold_search(model_func, pulse, Fs, min_I, max_I)
    max_iters = 25; % 25 iters guarantees extremely high precision
    stim = [pulse, zeros(1, round(5e-3 * Fs))]; 
    p_n = zeros(1, length(stim)); % Zero noise for deterministic search
    c_n = zeros(1, length(stim)); 
    
    for iter = 1:max_iters
        test_amp = (min_I + max_I) / 2;
        [~, SpT, ~, ~, ~, ~] = model_func(stim * test_amp, p_n, c_n, Fs);
        if ~isempty(SpT)
            max_I = test_amp; 
        else
            min_I = test_amp; 
        end
    end
    threshold = (min_I + max_I) / 2;
end

function latency_us = get_deterministic_latency(model_func, pulse, amp, Fs)
    stim = [pulse * amp, zeros(1, round(5e-3 * Fs))]; 
    p_n = zeros(1, length(stim)); 
    c_n = zeros(1, length(stim)); 
    [~, SpT, ~, ~, ~, ~] = model_func(stim, p_n, c_n, Fs);
    if ~isempty(SpT)
        latency_us = SpT(1) * 1e6;
    else
        latency_us = NaN;
    end
end