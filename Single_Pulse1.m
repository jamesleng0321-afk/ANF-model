clear; clc; close all;

%% Model Parameters %%
Fs = 1e6; % Sampling frequency 
NoiseAlpha = 0.8; 

%% Cathodic pulse (阴极脉冲 - 论文对应蓝色倒三角) %%
PPD = 39; % 论文 Fig 3 采用的是 39 微秒
Istim = [0, -ones(1,PPD), zeros(1,1000)]; 

[Level,Probability,Latency1,Jitter1]=Library.FindThreshold(Istim,Fs,NoiseAlpha,0.0001e-6,@Model_SinglePulse,1000);

% 数据单位转换：Amps -> mA, Seconds -> us
Level_mA = Level * 1000; 
Latency1_us = Latency1 * 1e6;
Jitter1_us = Jitter1 * 1e6;

figure('Position', [100, 100, 400, 700]);

% 1. 绘制放电概率 (FE Curve)
subplot(3,1,1)
plot(Level_mA, Probability, 'v', 'MarkerFaceColor', '#5B84C4', 'MarkerEdgeColor', '#5B84C4', 'MarkerSize', 8); hold on;
[muCathode,sigma,xtemp,ytemp]=Library.FitNeuronDynamicRange(Level',Probability);
plot(xtemp*1000, ytemp, '-', 'Color', '#5B84C4', 'LineWidth', 1.5);

% 2. 绘制潜伏期 (Latency)
subplot(3,1,2)
plot(Level_mA, Latency1_us, 'v', 'MarkerFaceColor', '#5B84C4', 'MarkerEdgeColor', '#5B84C4', 'MarkerSize', 8); hold on;

% 3. 绘制抖动 (Jitter)
subplot(3,1,3)
plot(Level_mA, Jitter1_us, 'v', 'MarkerFaceColor', '#5B84C4', 'MarkerEdgeColor', '#5B84C4', 'MarkerSize', 8); hold on;

%% Anodic pulse (阳极脉冲 - 论文对应绿色正三角) %%
pause(0.001);
[Level,Probability,Latency2,Jitter2]=Library.FindThreshold(-Istim,Fs,NoiseAlpha,0.0001e-6,@Model_SinglePulse,1000);

% 数据单位转换：Amps -> mA, Seconds -> us
Level_mA = Level * 1000; 
Latency2_us = Latency2 * 1e6;
Jitter2_us = Jitter2 * 1e6;

% 1. 补充放电概率
subplot(3,1,1)
plot(Level_mA, Probability, '^', 'MarkerFaceColor', '#74C374', 'MarkerEdgeColor', '#74C374', 'MarkerSize', 8); 
[muAnode,sigma,xtemp,ytemp]=Library.FitNeuronDynamicRange(Level',Probability);
plot(xtemp*1000, ytemp, '-', 'Color', '#74C374', 'LineWidth', 1.5);
xlim([0.4, 0.9]);
ylim([0, 1.05]);
ylabel('Probability of spiking');
box off; set(gca, 'TickDir', 'out');

% 2. 补充潜伏期
subplot(3,1,2)
plot(Level_mA, Latency2_us, '^', 'MarkerFaceColor', '#74C374', 'MarkerEdgeColor', '#74C374', 'MarkerSize', 8); 
xlim([0.4, 0.9]);
ylim([150, 600]);
ylabel('Latency (\mu s)');
box off; set(gca, 'TickDir', 'out');

% 3. 补充抖动
subplot(3,1,3)
plot(Level_mA, Jitter2_us, '^', 'MarkerFaceColor', '#74C374', 'MarkerEdgeColor', '#74C374', 'MarkerSize', 8);
xlim([0.4, 0.9]);
ylim([0, 200]);
ylabel('Jitter (\mu s)');
xlabel('Pulse level (mA)');
box off; set(gca, 'TickDir', 'out');

% 统一添加图例
subplot(3,1,1);
legend('Cathodic', '', 'Anodic', '', 'Location', 'southeast');