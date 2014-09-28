% Phase 1 ECE 4784
% Kyle Francis
clear all;
clc;
%% Creating timestep %%
ts = input('Enter a timestep in ms: ');
time = 100; %ms
x = 1:ts:time;
len = length(x);

%% Constants %%
g_k = 36; %mS/cm^2   maximum conductances
g_na = 120; %mS/cm^2
g_l = 0.3; %mS/cm^2
E_k = -12; %mV
E_na = 115; %mV
E_l = 10.6; %mV
V_rest = -70; %mV
Vm(1) = 0; %V
Cm = 1.0; %uF/cm^2

%% Gating variables %%
am(1) = 0.1*((25-Vm(1))/(exp((25-Vm(1))/10) - 1));
bm(1) = 4*exp(-Vm(1)/18);
an(1) = 0.01*((10-Vm(1))/(exp((10-Vm(1))/10) - 1));
bn(1) = 0.125*exp(-Vm(1)/80);
ah(1) = 0.07*exp(-Vm(1)/20);
bh(1) = 1/(exp((30-Vm(1))/10) + 1);

m(1) = am(1)/(am(1)+bm(1));
n(1) = an(1)/(an(1)+bn(1));
h(1) = ah(1)/(ah(1)+bh(1));

%% Currents %%
% Current = activation probability (m) [or inactivation prob (n)] * max conduct * Nernst potential
I_na(1) = m(1)^3*g_na*h(1)*(Vm(1)-E_na);
I_k(1) = n(1)^4*g_k*(Vm(1)-E_k);
I_l(1) = g_l*(Vm(1)-E_l);
amp = input('Enter the amplitude in mA of the injected current: ');
dur = input('Enter the duration of the injected current [between 0 and 100 ms]: ');
I = zeros(1, len);
for x = 1:dur
    I(x) = amp;
end
I_ion(1) = I(1) - I_k(1) - I_na(1) - I_l(1);

%% Derivatives %%
dVm = I_ion/Cm;
dm = am*(1-m)-bm*m;
dn = an*(1-n)-bn*n;
dh = ah*(1-h)-bh*h;

%% Euler's method %%
%yn+1 = yn+h*dy/dt;
for i = 1:len-1
    am(i) = 0.1*((25-Vm(i))/(exp((25-Vm(i))/10) - 1));
    bm(i) = 4*exp(-Vm(i)/18);
    an(i) = 0.01*((10-Vm(i))/(exp((10-Vm(i))/10) - 1));
    bn(i) = 0.125*exp(-Vm(i)/80);
    ah(i) = 0.07*exp(-Vm(i)/20);
    bh(i) = 1/(exp((30-Vm(i))/10) + 1);
    
    dm(i) = am(i)*(1-m(i))-(bm(i)*m(i));
    dn(i) = an(i)*(1-n(i))-bn(i)*n(i);
    dh(i) = ah(i)*(1-h(i))-bh(i)*h(i);
    
    m(i+1) = m(i) + ts*dm(i);
    n(i+1) = n(i) + ts*dn(i);
    h(i+1) = h(i) + ts*dh(i);
    
    I_na(i) = m(i).^3*g_na*h(i)*(Vm(i)-E_na);
    I_k(i) = n(i).^4*g_k*(Vm(i)-E_k);
    I_l(i) = g_l*(Vm(i)-E_l);
    I_ion(i) = I(i) - I_k(i) - I_na(i) - I_l(i);
    
    dVm(i) = I_ion(i)/Cm;
    Vm(i+1) = Vm(i) + ts*dVm(i);
    
end
%% Plotting %%
Vm = Vm + V_rest;
plot(x,Vm)
title('Phase 1');
xlabel('Time (ms)');
ylabel('Voltage potential (mV)');