%% ========================================================================
% Heffer2010_psBLIF_vs_ANF.m
% Reproduction of Heffer 2010 rate-dependent spike probability changes
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% Load experimental data
% ------------------------------------------------------------------------
load("Heffer2010.mat"); % Ensure this file exists in your path
stim_rates = [200, 1000, 2000, 5000];    % pulses per second
signal_length_s = 2e-3;                  % 2 ms (onset period)

%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------
phase_len_us = 25;
ipg_us       = 8;
Fs           = 1e6; % 1 MHz sampling rate for ANF

% ---- psBLIF Pulse ----
phase_len = phase_len_us * 1e-6;
ipg       = ipg_us * 1e-6;
pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = ipg;
pulse_psBLIF.negative_duration  = phase_len;
pulse_psBLIF.negative_amplitude = -1;

% ---- ANF Pulse ----
% [EN] Create single sample pulse for ANF (assuming 1us = 1 sample at 1MHz)
% [ZH] 创建 ANF 的单脉冲波形 (假设 1MHz 下 1微秒 = 1个采样点)
pulse_ANF = [0, -1 * ones(1, phase_len_us), ...
             zeros(1, ipg_us), ...
             1 * ones(1, phase_len_us), 0];

%% ------------------------------------------------------------------------
% Target probability levels / 目标概率区间
% ------------------------------------------------------------------------
low_probs    = (0.02:0.02:0.18);
medium_probs = (0.3:0.05:0.7);
high_probs   = (0.75:0.025:0.95);

%% ------------------------------------------------------------------------
% Compute single pulse thresholds (psBLIF + ANF)
% 计算单脉冲的电流幅度阈值
% ------------------------------------------------------------------------
% ================== psBLIF ==================
C = psBLIF_default_parameters();
single_ps(1) = pulse_psBLIF;
eval_single_ps = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_probe(a, single_ps, C), amp);

disp('Calculating psBLIF thresholds...');
low_levels_ps    = arrayfun(@(p) find_threshold(eval_single_ps, p), low_probs);
medium_levels_ps = arrayfun(@(p) find_threshold(eval_single_ps, p), medium_probs);
high_levels_ps   = arrayfun(@(p) find_threshold(eval_single_ps, p), high_probs);

% ================== ANF =====================
NoiseAlpha = 0.8;
disp('Calculating ANF dynamic range & thresholds...');
% [EN] Find threshold curve for a single pulse
% [ZH] 寻找单脉冲的阈值曲线
[Level_curve, Prob_curve] = Library.FindThreshold([pulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, sigmaSingle]   = Library.FitNeuronDynamicRange(Level_curve', Prob_curve);

% [EN] Use Inverse Gaussian CDF to find exact amplitudes for target probabilities
% [ZH] 使用高斯累积分布函数的逆函数，反推出目标概率对应的精确电流幅度
% Note: norminv(p, mu, sigma) = mu + sigma * sqrt(2) * erfinv(2*p - 1)
get_ANF_level = @(p) muSingle + sigmaSingle * sqrt(2) * erfinv(2*p - 1);

low_levels_ANF    = arrayfun(get_ANF_level, low_probs);
medium_levels_ANF = arrayfun(get_ANF_level, medium_probs);
high_levels_ANF   = arrayfun(get_ANF_level, high_probs);

%% ------------------------------------------------------------------------
% Preallocate / 预分配内存
% ------------------------------------------------------------------------
prob_change_ps.low    = nan(length(stim_rates), length(low_levels_ps));
prob_change_ps.medium = nan(length(stim_rates), length(medium_levels_ps));
prob_change_ps.high   = nan(length(stim_rates), length(high_levels_ps));

prob_change_ANF = prob_change_ps;

%% ------------------------------------------------------------------------
% Loop over stimulation rates / 遍历刺激频率
% ------------------------------------------------------------------------
nTrials_ANF = 50; % [EN] Monte Carlo trials. Increase to 100+ for smoother curves. / [ZH] 蒙特卡洛试验次数，增大可使曲线更平滑。

for s_ind = 1:length(stim_rates)
    rate = stim_rates(s_ind);
    fprintf('\nProcessing Rate: %d pps...\n', rate);
    
    % [EN] Calculate how many pulses fall into the 2ms onset window
    % [ZH] 计算在 2ms 的起始时间窗内有多少个脉冲
    ipi_s = 1 / rate;
    num_pulses_in_onset = 1 + floor(signal_length_s / ipi_s);
    
    % ---- Build psBLIF pulse train ----
    pulse_train_ps = build_rate_psBLIF_train(signal_length_s, rate, pulse_psBLIF);
    
    %% ---- LOW LEVELS ----
    prob_change_ps.low(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, low_levels_ps, C);
    prob_change_ANF.low(s_ind,:) = ...
        compute_rate_change_ANF(pulse_ANF, low_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_ANF, num_pulses_in_onset, low_probs);
        
    %% ---- MEDIUM LEVELS ----
    prob_change_ps.medium(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, medium_levels_ps, C);
    prob_change_ANF.medium(s_ind,:) = ...
        compute_rate_change_ANF(pulse_ANF, medium_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_ANF, num_pulses_in_onset, medium_probs);
        
    %% ---- HIGH LEVELS ----
    prob_change_ps.high(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, high_levels_ps, C);
    prob_change_ANF.high(s_ind,:) = ...
        compute_rate_change_ANF(pulse_ANF, high_levels_ANF, rate, signal_length_s, Fs, NoiseAlpha, nTrials_ANF, num_pulses_in_onset, high_probs);
end

%% ------------------------------------------------------------------------
% Plot / 绘图
% ------------------------------------------------------------------------
visual_off = 0.07;
fig = figure("Name","Heffer2010_comparison_ANF","DefaultAxesFontSize",13);
fig.OuterPosition(3:4) = [1500 600];
tiledlayout(1,3);

plot_panel(data_Heffer2010(3), prob_change_ps.high, ...
           prob_change_ANF.high, stim_rates, visual_off, "high");
plot_panel(data_Heffer2010(2), prob_change_ps.medium, ...
           prob_change_ANF.medium, stim_rates, visual_off, "medium");
plot_panel(data_Heffer2010(1), prob_change_ps.low, ...
           prob_change_ANF.low, stim_rates, visual_off, "low");

lg = legend("data", "psBLIF", "ANF", 'numcolumns', 3);
lg.Layout.Tile = "North";

%% ========================================================================
% Helper Functions / 辅助函数
% ========================================================================

% [EN] Wrapper to compute facilitation for stochastic ANF model
% [ZH] 计算随机 ANF 模型易化作用(Facilitation)的封装函数
function change = compute_rate_change_ANF(pulse, amplitude_levels, rate, signal_len, Fs, NoiseAlpha, nTrials, num_pulses, target_probs)
    change = zeros(1, length(amplitude_levels));
    
    % [EN] Generate the normalized pulse train sequence
    % [ZH] 生成归一化的脉冲序列
    Istim = Experiment.stim_PulseTrain(pulse, rate, 100, 0, signal_len, Fs);
    
    for i = 1:length(amplitude_levels)
        amp = amplitude_levels(i);
        input = Istim * amp;
        
        spikes_in_onset = 0;
        
        % [EN] Monte Carlo Trials
        % [ZH] 蒙特卡洛随机试验
        for tr = 1:nTrials
            p_noise = Library.oneonfnoise(length(input), NoiseAlpha);
            c_noise = Library.oneonfnoise(length(input), NoiseAlpha);
            [~, SpTimes, ~, ~] = Model_PulseTrain(input, p_noise, c_noise, Fs);
            
            % [EN] If at least one spike occurred within the onset window
            % [ZH] 如果在 onset 时间窗 (2ms) 内发生了至少一次脉冲
            if any(SpTimes <= signal_len)
                spikes_in_onset = spikes_in_onset + 1;
            end
        end
        
        % [EN] Measured probability is the fraction of trials with a spike
        % [ZH] 测量到的概率是产生脉冲的试验比例
        p_onset_measured = spikes_in_onset / nTrials;
        
        % [EN] Predicted probability based on single pulse (No facilitation)
        % [ZH] 基于单脉冲的预测概率 (假设无易化作用)
        p_first = target_probs(i);
        predicted = min(1, p_first * num_pulses);
        
        % [EN] Facilitation = Measured - Predicted
        % [ZH] 易化作用 = 实测值 - 预测值
        change(i) = p_onset_measured - predicted;
    end
end

function pulse_train = build_rate_psBLIF_train(signal_len, rate, pulse)
    ipi = 1 / rate;
    n   = 1 + floor(signal_len / ipi);
    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1) * ipi;
    end
end

function thr = find_threshold(eval_function, target_probability)
    thr = fzero(@(a) eval_function(a)-target_probability, [0,1]);
end

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
end

function [history, alifp_ret] = psBLIF_wrapper_scale_probe(amplitude, pulse_train, C)
    scaled = pulse_train;
    last = numel(scaled);
    scaled(last).positive_amplitude = scaled(last).positive_amplitude * amplitude;
    scaled(last).negative_amplitude = scaled(last).negative_amplitude * amplitude;
    [history, alifp_ret] = psBLIF(scaled, C);
end

function [history, alifp_ret] = psBLIF_wrapper_scale_all(amplitude, pulse_train, C)
    scaled = pulse_train;
    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = scaled(k).positive_amplitude * amplitude;
        scaled(k).negative_amplitude = scaled(k).negative_amplitude * amplitude;
    end
    [history, alifp_ret] = psBLIF(scaled, C);
end

function change = compute_rate_change_psBLIF(pulse_train, levels, C)
    change = zeros(1,length(levels));
    for i = 1:length(levels)
        level = levels(i);
        [~, ps_out] = psBLIF_wrapper_scale_all(level, pulse_train, C);
        per_pulse_prob = psBLIF_per_pulse_prob(ps_out);
        p_first = per_pulse_prob(1); 
        p_total = 1-prod(1-per_pulse_prob);
        change(i) = p_total - min(1, p_first * numel(pulse_train));
    end
end

% [EN] Helper function to calculate per-pulse prob for psBLIF
% [ZH] 辅助提取 psBLIF 每次脉冲概率的函数
function probs = psBLIF_per_pulse_prob(alifp_ret)
    probs = zeros(1, length(alifp_ret));
    for i = 1:length(alifp_ret)
        probs(i) = sum([alifp_ret{i}.path_prob]);
    end
end

function plot_panel(data_struct, ps_change, a_change, rates, off, label)
    nexttile
    lower_error = data_struct.change(1:4,2) - data_struct.change(5:2:end,2);
    upper_error = data_struct.change(6:2:end,2) - data_struct.change(1:4,2);
    errorbar((1:4)-off, data_struct.change(1:4,2), ...
             lower_error, upper_error, ...
             'square','linewidth',2,"color","black");
    hold on
    median_ps = median(ps_change,2);
    q25_ps = prctile(ps_change,25,2);
    q75_ps = prctile(ps_change,75,2);
    errorbar((1:4)+off, median_ps, ...
             median_ps-q25_ps, q75_ps-median_ps, ...
             'square','linewidth',2,"color","#0072BD");
             
    median_a = median(a_change,2);
    q25_a = prctile(a_change,25,2);
    q75_a = prctile(a_change,75,2);
    errorbar((1:4)+2.5*off, median_a, ...
            median_a-q25_a, q75_a-median_a, ...
            'square','linewidth',2,"color","#D95319");
            
    ylabel("Spike probability change")
    xticks(1:4)
    xticklabels(string(rates))
    xlabel("Pulse rate [pps]")
    title(label)
    ylim([-0.4 1])
    xlim([0.5 4.5])
    yline(0)
end