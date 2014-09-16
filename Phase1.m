% Phase 1 ECE 4784
% Kyle Francis


%Constants
time = 100e-3;
g_k = 36e-3; %S/cm^2
g_na = 120e-3; %S/cm^2
g_cl = 0.3e-3; %S/cm^2
E_k = -12e-3; %V
E_na = 115e-3; %V
E_cl = 10.6e-3; %V
V_rest = -70e-3; %V

%Gating variables
am = 0.1*((25-Vm)/(exp((25-Vm)/10) - 1));
bm = 4*exp(-Vm/18);
an = 0.01*((10-Vm)/(exp((10-Vm)/10) - 1));
bn = 0.125*exp(-Vm/80);
ah = 0.07*exp(-Vm/20);
bh = 1/(exp((30-Vm)/10) + 1);

%Currents
I_na = m^3*g_na*h*(Vm-E_na);
I_k = n^4*g_k*(Vm-E_k);
I_cl = g_cl*(Vm-E_cl);
I_ion = I - I_k - I_na - I_cl;

%Derivatives
dVm = I_ion/Cm;
dm = am*(1-m)-bm*m;
dn = an*(1-n)-bn*n;
dh = ah*(1-h)-bh*h;

%Euler's method
yn+1 = yn+h*f(tn,yn);
