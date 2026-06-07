function [Spike, SpTimes, p_sptimes, c_sptimes, p_v, c_v] = Model_PulseTrain_Bruce(Istim, p_Noise, c_Noise, Fs)
%% Two point neuron model function (Improved Parameterization)
% Based on Joshi et al. (2017) and modified according to Bruce (2024):
% "Improved Parameterization of a Hybrid Phenomenological-Biophysical Model..."

%% Parameters %%
%%%%%%%%%%%%%%%%%%% Common parameters %%%%%%%%%%%%%%%%%%%
El = -0.0800;
Vt = -0.0700;
Vr = -0.0840;
vPeak = 0.024;
AbsRef = 5.0000e-04;
a1 = 0.0026;
a2 = 0.005;
b = 90e-06;
InhibitAlpha =  0.75;

%%%%%%%%%%%%%%%%%%% Peripheral node %%%%%%%%%%%%%%%%%%%
p_Cm = 856.96e-09;
p_gL = 0.0011;
p_Dt = 0.010;
p_tauw1 = 400e-06;
p_tauw2 = 4500e-06;
% [Modification 1]: According to Bruce (2024), scale down RS by a factor of 3 
% to reduce the spontaneous firing rate to physiological levels.
p_RS = 0.062 / 3; 
p_Threshold = 543e-6;
p_Sigma = p_RS * p_Threshold;

%%%%%%%%%%%%%%%%%%% Central node %%%%%%%%%%%%%%%%%%%
c_Cm = 1772.4e-09;
c_gL = 0.0027;
c_Dt = 0.0030;
c_tauw1 = 250e-06;
c_tauw2 = 3000e-06;
% [Modification 1]: According to Bruce (2024), scale down RS by a factor of 3.
c_RS = 0.075 / 3; 
c_Threshold = 731e-6;
c_Sigma = c_RS * c_Threshold;

%% Model Computations %%
dt = 1/Fs;
Spike = 0;
SpTimes = [];
nt = length(Istim);

% [Modification 2]: Randomization of membrane potential initialization 
% (Reference: Bruce 2024 Fig 4)
% Original standard deviations were 4.77 mV and 3.23 mV. 
% The standard deviation scales linearly with the RS parameter.
std_p = 4.77e-3 / 3; % Scaled peripheral std dev is approx 1.59 mV
std_c = 3.23e-3 / 3; % Scaled central std dev is approx 1.08 mV

% Peripheral node setup
p_I = Istim * -1 ;
p_I(p_I < 0) = InhibitAlpha * p_I(p_I < 0);
p_t = 0:dt:length(p_I)/Fs;
p_v = zeros(1, length(p_t)); 
p_v(:) = El;  
% Introduce normally-distributed initial random offset
p_v(1) = El + std_p * randn(1); 

p_sptimes = [];
p_Noise = p_Sigma * p_Noise;

% Central node setup
c_I = Istim;
c_I(c_I < 0) = InhibitAlpha * c_I(c_I < 0);
c_t = 0:dt:length(c_I)/Fs;
c_v = zeros(1, length(c_t)); 
c_v(:) = El;  
% Introduce normally-distributed initial random offset
c_v(1) = El + std_c * randn(1); 

c_sptimes = [];
c_Noise = c_Sigma * c_Noise;

TimeSinceSpike = 1000;
SpikingNow = 0;
ip = 0;

% Initialize adaptive currents
p_w1 = zeros(1, length(p_t)); 
p_w2 = zeros(1, length(p_t));
c_w1 = zeros(1, length(c_t)); 
c_w2 = zeros(1, length(p_t)); 

for idt = 2:nt
   
    % Calculate adaptation-triggered adaptive current
    p_w1(idt) = p_w1(idt-1) + dt/p_tauw1*(a1*(p_v(idt-1)-El) - p_w1(idt-1));
    p_w2(idt) = p_w2(idt-1) + dt/p_tauw2*(a2*(p_v(idt-1)-El) - p_w2(idt-1));
    c_w1(idt) = c_w1(idt-1) + dt/c_tauw1*(a1*(c_v(idt-1)-El) - c_w1(idt-1));
    c_w2(idt) = c_w2(idt-1) + dt/c_tauw2*(a2*(c_v(idt-1)-El) - c_w2(idt-1));
    
    TimeSinceSpike = TimeSinceSpike + dt;
    if TimeSinceSpike < AbsRef %% If the neuron is in the absolute refractory period
        p_v(idt) = (p_v(idt-1) + ...
            dt/p_Cm*(p_gL*(El-p_v(idt-1))  + ...
            p_gL*p_Dt*exp((p_v(idt-1)-Vt)/p_Dt) - ...
            p_w1(idt-1) - p_w2(idt-1) + p_Noise(idt)));
        
        c_v(idt) = (c_v(idt-1) + ...
            dt/c_Cm*(c_gL*(El-c_v(idt-1))  + ...
            c_gL*c_Dt*exp((c_v(idt-1)-Vt)/c_Dt) - ...
            c_w1(idt-1) - c_w2(idt-1) + c_Noise(idt)));
        
    else %% If neuron is not in the absolute refractory period
        p_v(idt) = (p_v(idt-1) + ...
            dt/p_Cm*(p_gL*(El-p_v(idt-1))  + ...
            p_gL*p_Dt*exp((p_v(idt-1)-Vt)/p_Dt) - ...
            p_w1(idt-1) - p_w2(idt-1) + ...
            p_I(idt-1) + p_Noise(idt)));
        
        c_v(idt) = (c_v(idt-1) + ...
            dt/c_Cm*(c_gL*(El-c_v(idt-1))  + ...
            c_gL*c_Dt*exp((c_v(idt-1)-Vt)/c_Dt) - ...
            c_w1(idt-1) - c_w2(idt-1) + ...
            c_I(idt-1) + c_Noise(idt)));
        
        if p_v(idt) > vPeak % if spike occurs at peripheral node
            SpikingNow = 1;
            ip = idt-1 + (vPeak-p_v(idt-1))/(p_v(idt)-p_v(idt-1));   % Estimate spike time via interpolation
            p_sptimes = [p_sptimes, ip*dt];
        elseif c_v(idt) > vPeak % if spike occurs at central node
            SpikingNow = 1;
            ip = idt-1 + (vPeak-c_v(idt-1))/(c_v(idt)-c_v(idt-1));   % Estimate spike time via interpolation
            c_sptimes = [c_sptimes, ip*dt];
        else
            SpikingNow = 0;
        end
        
        if SpikingNow % Process spike event for the dual-node setup
            TimeSinceSpike = 0;
            SpTimes = [SpTimes, ip*dt];
            Spike = Spike + 1;
            
            % Reset the membrane voltage for both nodes
            p_v(idt) = Vr;
            c_v(idt) = Vr;
            
            % Update adaptation parameters for both nodes
            p_w1(idt) = p_w1(idt-1) + (ip-idt+1)*dt/p_tauw1*(a1*(p_v(idt-1)-El) - p_w1(idt-1));
            p_w2(idt) = p_w2(idt-1) + (ip-idt+1)*dt/p_tauw2*(a2*(p_v(idt-1)-El) - p_w2(idt-1)) + b;            
            c_w1(idt) = c_w1(idt-1) + (ip-idt+1)*dt/c_tauw1*(a1*(c_v(idt-1)-El) - c_w1(idt-1));
            c_w2(idt) = c_w2(idt-1) + (ip-idt+1)*dt/c_tauw2*(a2*(c_v(idt-1)-El) - c_w2(idt-1)) + b;

            SpikingNow = 0;
        end
    end
end
return