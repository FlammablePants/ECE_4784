% Phase 1 ECE 4784
% Kyle Francis
clear all;
clc;

%% Creating timestep %%
ts = 0.05; % if too large, the derivatives in the Eulers formula would be too large and cause an exp growth
time = 100; %ms
x = 0 : ts : time;
len = length(x); % time array length

%% Constants %%
g_kbar = 36; %mS/cm^2   maximum conductances
g_nabar = 120; %mS/cm^2
g_lbar = 0.3; %mS/cm^2
E_k = -12; %mV
E_na = 115; %mV
E_l = 10.6; %mV
V_rest = -70; %mV
Cm = 1.0; %uF/cm^2

%% Gating variables %%
Vm(1) = 0; %V
am(1) = 0.1 * ((25 - Vm(1)) / (exp((25 - Vm(1)) / 10) - 1));
bm(1) = 4 * exp(-Vm(1) / 18);
an(1) = 0.01 * ((10 - Vm(1)) / (exp((10 - Vm(1)) / 10) - 1));
bn(1) = 0.125 * exp(-Vm(1) / 80);
ah(1) = 0.07 * exp(-Vm(1) / 20);
bh(1) = 1 / (exp((30 - Vm(1)) / 10) + 1);

m(1) = am(1) / (am(1) + bm(1));
n(1) = an(1) / (an(1) + bn(1));
h(1) = ah(1) / (ah(1) + bh(1));

%% Injected Current %%
amp = input('Enter the amplitude in mA of the injected current: ');
dur = input('Enter the duration of the injected current [between 0 and 100 ms]: ');
I = zeros(1, len); % Setting up I to be the same size as the time array
for y = 1:floor(dur/ts) % injected current goes from 0 to (duration entered / time step), elsewhere is 0
    I(y) = amp;
end

%% Euler's method %%
%yn+1 = yn+h*dy/dt;

for i = 1 : len-1
    am = 0.1 * ((25 - Vm(i)) / (exp((25 - Vm(i)) / 10) - 1));
    bm = 4 * exp(-Vm(i) / 18);
    an = 0.01 * ((10 - Vm(i)) / (exp((10 - Vm(i)) / 10) - 1));
    bn = 0.125 * exp(-Vm(i) / 80);
    ah = 0.07 * exp(-Vm(i) / 20);
    bh = 1 / (exp((30 - Vm(i)) / 10) + 1);
    
    g_na(i) = m(i)^3 * g_nabar * h(i);
    g_k(i) = n(i)^4 * g_kbar;
    g_l(i) = g_lbar;
    
    I_na(i) = g_na(i) * (Vm(i)-E_na);
    I_k(i) = g_k(i) * (Vm(i)-E_k);
    I_l(i) = g_l(i) * (Vm(i)-E_l);
    I_ion(i) = I(i) - I_k(i) - I_na(i) - I_l(i);
    
    dm = am * (1-m(i)) - bm * m(i);
    dn = an * (1-n(i)) - bn * n(i);
    dh = ah * (1-h(i)) - bh * h(i);
    
    m(i+1) = m(i) + ts * dm;
    n(i+1) = n(i) + ts * dn;
    h(i+1) = h(i) + ts * dh;
    
    dVm(i) = I_ion(i) / Cm;
    Vm(i+1) = Vm(i) + ts * dVm(i);
    
end

%% Plotting %%
Vm = Vm + V_rest; % setting Vm to be relative to V_rest
figure(1)
plot(x,Vm)
title('Membrane Potential');
xlabel('Time (ms)');
ylabel('Membrane Voltage (mV)');
legend('Vm');
axis([0, 100, -100, 80]);

%Plot conductances
newtime = x(1:len-1); % Creating time to be same length as conductances
figure(2)
plot(newtime, g_na, 'r', newtime, g_k, 'b')
title('gK and gNa')
xlabel('Time (ms)')
ylabel('Conductance (mS/cm^2)')
legend('gNa','gK')

%Plot