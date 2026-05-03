%% THROTTLE VALVE — JOULE-THOMSON TEMPERATURE DROP
% Reactor outlet gas throttled from 200 bar to 80 bar
% Purpose: bring pressure BELOW NH3 critical pressure (Pc = 113.5 bar)
% so that NH3 behaves as a normal subcritical fluid in the condenser.
% Throttling is ISENTHALPIC — no heat exchange duty, no work.
% Temperature changes due to the Joule-Thomson effect.
%
% Joule-Thomson coefficient: mu_JT = (dT/dP) at constant H
% For real gases: mu_JT = (1/Cp) * [T*(dV/dT)_P - V]
% For ideal gas: mu_JT = 0 (no temperature change)
% N2 and H2 are nearly ideal at these conditions → small mu_JT
% NH3 is more non-ideal but at high T it also approaches ideal behaviour
% We use mixture mu_JT estimated from literature values at these conditions.

clc; clear;

%% ============================================================
%  WHY 200 bar -> 80 bar?
%  NH3 critical pressure Pc = 113.5 bar
%  NH3 critical temperature Tc = 132.4 deg C (405.5 K)
%  At 200 bar (> Pc), NH3 is supercritical — no defined boiling point,
%  no true latent heat, Antoine equation invalid.
%  Throttling to 80 bar (< Pc) puts NH3 below its critical pressure,
%  restoring normal vapour-liquid behaviour for the condenser.
%  80 bar chosen as a round number comfortably below 113.5 bar,
%  high enough to keep good NH3 condensation temperature (avoids
%  needing very deep refrigeration) and manageable compressor ratio.
%% ============================================================

%% Inlet conditions (reactor outlet)
T1 = 450 + 273.15;   % K  — reactor operating temperature
P1 = 200;            % bar — reactor operating pressure
P2 = 80;             % bar — target: below NH3 Pc = 113.5 bar
dP = P2 - P1;        % bar — negative (pressure drop)

%% Molar flow rates (kmol/day)
F_NH3 = 3095.98;
F_N2  = 5882.35;
F_H2  = 17647.06;
F_total = F_NH3 + F_N2 + F_H2;

%% Mole fractions
y_NH3 = F_NH3 / F_total;
y_N2  = F_N2  / F_total;
y_H2  = F_H2  / F_total;

%% Joule-Thomson coefficients at inlet conditions (~450 deg C, 200 bar)
% Units: K/bar
% N2:  near ideal at high T, mu_JT ~ +0.01 K/bar (very small, slightly positive)
% H2:  inversion temperature ~200 K, so at 450 deg C mu_JT is NEGATIVE
%      meaning H2 warms slightly on expansion — mu_JT ~ -0.03 K/bar
% NH3: more non-ideal, mu_JT ~ +0.15 K/bar at these conditions
%      (still positive — cools on expansion)
% Sources: Smith Van Ness Abbott Appendix B; Perry's Chemical Engineers Handbook

mu_JT_NH3 =  0.15;    % K/bar
mu_JT_N2  =  0.01;    % K/bar
mu_JT_H2  = -0.03;    % K/bar  (negative: H2 heats on expansion above inversion T)

%% Mixture Joule-Thomson coefficient (mole fraction weighted)
mu_JT_mix = y_NH3*mu_JT_NH3 + y_N2*mu_JT_N2 + y_H2*mu_JT_H2;

%% Outlet temperature after throttling
% T2 = T1 + mu_JT_mix * dP   (dP is negative)
T2 = T1 + mu_JT_mix * dP;   % K
T2_C = T2 - 273.15;          % deg C

%% Throttle duty (always zero — isenthalpic process)
Q_throttle = 0;   % kW

%% ============================================================
%  DISPLAY RESULTS
%% ============================================================

fprintf('=======================================================\n');
fprintf('           THROTTLE VALVE RESULTS                     \n');
fprintf('=======================================================\n');
fprintf('NH3 critical pressure Pc   : 113.5 bar\n');
fprintf('NH3 critical temperature Tc: 132.4 deg C\n');
fprintf('-------------------------------------------------------\n');
fprintf('Inlet pressure  (P1)       : %7.1f bar\n', P1);
fprintf('Outlet pressure (P2)       : %7.1f bar  (<< Pc)\n', P2);
fprintf('Inlet temperature  (T1)    : %7.2f K  (%6.2f deg C)\n', T1, T1-273.15);
fprintf('-------------------------------------------------------\n');
fprintf('Joule-Thomson coeff - NH3  : %7.4f K/bar\n', mu_JT_NH3);
fprintf('Joule-Thomson coeff - N2   : %7.4f K/bar\n', mu_JT_N2);
fprintf('Joule-Thomson coeff - H2   : %7.4f K/bar  (negative: heats on expansion)\n', mu_JT_H2);
fprintf('Mixture mu_JT              : %7.4f K/bar\n', mu_JT_mix);
fprintf('-------------------------------------------------------\n');
fprintf('Outlet temperature (T2)    : %7.2f K  (%6.2f deg C)\n', T2, T2_C);
fprintf('Temperature drop           : %7.2f K\n', T2-T1);
fprintf('-------------------------------------------------------\n');
fprintf('Throttle heat duty         : %7.1f kW  (isenthalpic — no duty)\n', Q_throttle);
fprintf('=======================================================\n');
fprintf('\nFeed to condenser: %.2f deg C, %.0f bar\n', T2_C, P2);