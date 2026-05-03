%% =========================================================
%  ADIABATIC THROTTLING - ENERGY BALANCE
%  Process: H2 product from PSA
%  Inlet:   30 bar, 27 deg C
%  Outlet:  10 bar, T_out = ?
%  Assumption: Adiabatic, steady-state, no work
%  h1 = h2 (isenthalpic process)
%% =========================================================

clc; clear; close all;

%% ---- INLET CONDITIONS ----------------------------------
P1 = 30;           % bar
T1 = 27 + 273.15;  % K (300.15 K)
P2 = 10;           % bar

%% ---- GAS PROPERTIES FOR H2 ----------------------------
% Cp of H2 at ~300 K (approx constant over small range)
Cp = 28.82 %J/(mol.K) 

% Joule-Thomson coefficient for H2 at ~300 K, ~20-30 bar
% mu_JT = (dT/dP)_h  [K/bar]
% H2 at room temp has NEGATIVE mu_JT => warms on expansion
% Typical value: -0.03 K/bar at 300 K
mu_JT = -0.03;     % K/bar  (negative = warming on pressure drop)

%% ---- JOULE-THOMSON RELATION ----------------------------
% dT = mu_JT * dP
% For finite change: T2 = T1 + mu_JT * (P2 - P1)

dP    = P2 - P1;           % bar  (negative = pressure drop)
dT    = mu_JT * dP;        % K
T2    = T1 + dT;           % K
T2_C  = T2 - 273.15;       % deg C

%% ---- ENTHALPY CHECK (should be zero for isenthalpic) ---
% dh = Cp * dT - mu_JT * Cp * dP   [for real gas: dh = Cp*(1 - mu_JT*Cp)*dP ... 
% For isenthalpic: dh = 0 by definition
% We verify using: dh = Cp * dT + (v - T*dv/dT)*dP
% Simplified real-gas check using JT definition:
% h2 - h1 = Cp*(T2-T1) - mu_JT*Cp*(P2-P1)  --> should = 0
dh = Cp * dT - mu_JT * Cp * dP;   % J/kg  (should be ~0)

%% ---- ENERGY BALANCE DISPLAY ----------------------------
fprintf('============================================\n');
fprintf('   ADIABATIC THROTTLING - ENERGY BALANCE   \n');
fprintf('============================================\n');
fprintf('INLET  CONDITIONS:\n');
fprintf('   Temperature  T1 = %.2f K  (%.2f deg C)\n', T1, T1-273.15);
fprintf('   Pressure     P1 = %.1f bar\n', P1);
fprintf('\nOUTLET CONDITIONS:\n');
fprintf('   Pressure     P2 = %.1f bar\n', P2);
fprintf('   Temperature  T2 = %.4f K  (%.4f deg C)\n', T2, T2_C);
fprintf('\nCHANGES:\n');
fprintf('   dP = %.1f bar\n', dP);
fprintf('   dT = %.4f K  (%.4f deg C)\n', dT, dT);
fprintf('   dh = %.4f J/kg  (should be ~0 for isenthalpic)\n', dh);
fprintf('\nJoule-Thomson coefficient:\n');
fprintf('   mu_JT = %.4f K/bar\n', mu_JT);
fprintf('   => H2 WARMS slightly on expansion (negative mu_JT)\n');
fprintf('\nEXTRA COOLING NEEDED:\n');
fprintf('   Target outlet temp = 27 deg C = 300.15 K\n');
fprintf('   T2 after valve     = %.4f deg C\n', T2_C);
fprintf('   Trim cooler duty   = %.4f deg C to remove\n', T2_C - 27);
fprintf('============================================\n');

%% ---- COOLER HEAT DUTY (per kg of H2) --------------
% Assume mass flow rate input
n_dot = 4411.7646;          % kmol/day  -- change to your actual flow rate
T_target = 300;  % K  (27 deg C)

Q_cooler = n_dot * Cp * (T2 - T_target);   % kmol/day J/mol.K K = kJ/day

fprintf('\nCOOLER DUTY ');
fprintf('   Q_cooler = %.4f kJ/day = %.4f MW\n', Q_cooler, Q_cooler/1000/86400);

fprintf('============================================\n');

