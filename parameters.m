function pars = parameters

%% Reversal potentials (mV)
pars.vL = -62.5;
pars.vNa = 45;
pars.vK = -105;
pars.vH = -35;
pars.vCa = 120;

%% Conductances (mS/cm^2)
% gL = 2.5;
% gNa = 29.17;
% gK = 12.96;
% gH = 20;
% gLVA = 11.0213;
% gNaP = 8.3244;
% gHVA = 2;
% gBK = 5;
% gNew = 5;

pars.gL = 2.5;
pars.gNa = 29.17;
pars.gK = 12.96;
pars.gH = 20;
pars.gLVA = 15.0213;
pars.gNaP = 8.3244;
pars.gHVA = 2.0;
pars.gBK = 5;
pars.gNew = 10;

%% Na
pars.theta_mNa = -25;
pars.sigma_mNa = -6.5;

%% K
pars.theta_nK = -26;
pars.sigma_nK = -9;
pars.tau_nK = 10;

%% LVA
pars.theta_mLVA = -37.1;
pars.sigma_mLVA = -4.8916;
pars.tau_mLVA = 40;
pars.theta_hLVA = -59.2;
pars.sigma_hLVA = 11.2326;
pars.tau_hLVA = 350;

% HVA
pars.theta_mHVA = -10.0;
pars.sigma_mHVA = -6.5;

%% NaP
pars.theta_mNaP = -40;
pars.sigma_mNaP = -4;
pars.theta_hNaP = -54;
pars.sigma_hNaP = 5;
pars.tau_hNaP = 500;

%% H
pars.theta_hH = -61.32;
pars.sigma_hH = 5.855;
pars.tau_hH_T = 100;
pars.delta_hH_T = 0.205;
pars.theta_hH_T = -65.95;
pars.sigma_hH_T = 4.44;

% New
pars.theta_mNew = -40.0;
pars.sigma_mNew = -2.0;
pars.theta_nNew = -30.0;
pars.sigma_nNew = -2.0;

%% BK
pars.wBK_base = 170; %base time constant (ms)
pars.sigma_wBK = -15.6;

%% other parameters
pars.Ca0 = 0.00002; %concentration of calcium (mM)
pars.tau_Ca = 8; %calcium concentration time constant (ms)
pars.Ca_buffer = 0.5; %accounts for quick calcium buffering
pars.Ca_z = 2; %unitless number
pars.d = 1; %depth where calcuim concentration is relevent (microns)
pars.C = 21; % capacitance (microF/cm^2)
pars.F = 96485; %Faraday's constant (C/mol)

end
