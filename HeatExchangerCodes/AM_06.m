%% CONDENSER DUTY CALCULATION
% Inlet : throttle valve outlet (~450 deg C, 80 bar)
% Outlet: T_out back-calculated from Antoine equation (see below)
% Components: NH3, N2, H2
% Cp/R = A + B*T + C*T^2 + D/T^2  (T in Kelvin)
%
% WHY 80 bar?
%   NH3 critical pressure Pc = 113.5 bar. Reactor runs at 200 bar
%   which is supercritical for NH3 — no defined boiling point there.
%   Throttle valve upstream drops to 80 bar (< Pc) so NH3 condenses normally.
%
% HOW IS T_out DETERMINED?
%   We want exactly 2.941 kmol/day NH3 to remain as vapor in the recycle gas
%   and 2941.18 kmol/day to condense as liquid product.
%   This fixes the NH3 mole fraction in the gas phase, which fixes the
%   NH3 partial pressure, which via Antoine equation fixes T_out.
%
%   Steps:
%   1. Total gas at separator top = N2 + H2 + NH3 vapor
%      = 5882.35 + 17647.06 + 2.941 = 23532.35 kmol/day
%   2. y_NH3 = 2.941 / 23532.35 = 1.2497e-4
%   3. P_NH3 = y_NH3 * P_total = 1.2497e-4 * 80 = 0.009998 bar = 7.499 mmHg
%   4. Antoine: T = B / (A - log10(P_mmHg)) - C
%      T = 926.13 / (6.8929 - log10(7.499)) - 240.17 = -86.3 deg C
%
%   

clc; clear;

%% ============================================================
%  USER INPUTS 
%% ============================================================

% Latent heat of vaporisation of NH3 at 80 bar (kJ/kmol)
% At 80 bar NH3 is subcritical — real latent heat exists.
% Saturation temp at 80 bar ~ 50 deg C.
% Typical value ~ 17,000-19,000 kJ/kmol at 80 bar.
% Enter your value from NIST / textbook appendix:
delta_H_vap = 18000;      % kJ/kmol  <-- ENTER YOUR VALUE at 80 bar

% Cp/R = A + B*T + C*T^2 + D/T^2   (T in Kelvin)
% R = 8.314 kJ/kmol.K
% Default: Smith Van Ness Abbott Appendix C

%        A         B           C           D
ABCD_NH3 = [3.578,  3.020e-3,   0.000,     -0.186e5];
ABCD_N2  = [3.280,  0.593e-3,   0.000,      0.040e5];
ABCD_H2  = [3.249,  0.422e-3,   0.000,      0.083e5];

%% ============================================================
%  FIXED PROCESS VALUES
%% ============================================================

R = 8.314;   % kJ/kmol.K

% Molar flow rates (kmol/day)
F_NH3_total = 3095.98;
F_N2        = 5882.35;
F_H2        = 17647.06;

% Convert to kmol/s
F_NH3_s = F_NH3_total / 86400;
F_N2_s  = F_N2        / 86400;
F_H2_s  = F_H2        / 86400;

% NH3 split specification
F_NH3_condensed = 2941.18;   % kmol/day — liquid product (bottom of separator)
F_NH3_vapor     = 2.941;     % kmol/day — stays as vapor in recycle gas (top)

% Verify split adds up
if abs((F_NH3_condensed + F_NH3_vapor) - F_NH3_total) > 0.1
    warning('NH3 split does not sum to total NH3 flow. Check values.');
end

% Condenser inlet temperature (throttle valve outlet)
% JT cooling from 200->80 bar at 450 deg C is negligible (~few deg C)
% so inlet is essentially 450 deg C.

T_in   = 450 + 273.15;   % K
P_cond = 80;             % bar

%% ============================================================
%  BACK-CALCULATE T_out FROM ANTOINE EQUATION
%
%  NH3 vapor in recycle = 2.941 kmol/day
%  This sets the NH3 partial pressure in the gas phase.
%  Antoine equation then gives T at which NH3 vapor pressure = P_NH3.
%
%  Antoine (P in mmHg, T in deg C):
%    log10(P) = A - B/(C + T)
%  Rearranged:
%    T = B / (A - log10(P_mmHg)) - C
%
%  Valid range: approx -33 to +50 deg C for this Antoine set.
%  Our result will be well below -33 deg C — Antoine is being
%  extrapolated here. 
%% ============================================================

F_gas_top  = F_N2 + F_H2 + F_NH3_vapor;        % total gas at separator top (kmol/day)
y_NH3_gas  = F_NH3_vapor / F_gas_top;           % NH3 mole fraction in gas phase
P_NH3_bar  = y_NH3_gas * P_cond;               % NH3 partial pressure (bar)
P_NH3_mmHg = P_NH3_bar * 750.062;              % convert to mmHg

% Antoine constants
Ant_A = 6.8929;
Ant_B = 926.13;
Ant_C = 240.17;

T_out_C = Ant_B / (Ant_A - log10(P_NH3_mmHg)) - Ant_C;   % deg C
T_out   = T_out_C + 273.15;                               % K

%% ============================================================
%  Cp INTEGRATION FUNCTION
%  Integral of (Cp/R) dT from Ta to Tb:
%  = R*[A*(Tb-Ta) + B/2*(Tb^2-Ta^2) + C/3*(Tb^3-Ta^3) - D*(1/Tb-1/Ta)]
%% ============================================================

function H = enthalpy_change(ABCD, Ta, Tb, R)
    A = ABCD(1); B = ABCD(2); C = ABCD(3); D = ABCD(4);
    H = R * ( A*(Tb - Ta) ...
            + (B/2)*(Tb^2 - Ta^2) ...
            + (C/3)*(Tb^3 - Ta^3) ...
            - D*(1/Tb - 1/Ta) );
    % kJ/kmol — negative when cooling (Tb < Ta)
end

%% ============================================================
%  SENSIBLE HEAT — cool all gas from T_in to T_out
%% ============================================================

dH_NH3 = enthalpy_change(ABCD_NH3, T_in, T_out, R);   % kJ/kmol
dH_N2  = enthalpy_change(ABCD_N2,  T_in, T_out, R);
dH_H2  = enthalpy_change(ABCD_H2,  T_in, T_out, R);

Q_sensible_NH3 = F_NH3_s * dH_NH3;   % kW
Q_sensible_N2  = F_N2_s  * dH_N2;
Q_sensible_H2  = F_H2_s  * dH_H2;

Q_sensible_total = Q_sensible_NH3 + Q_sensible_N2 + Q_sensible_H2;   % kW

%% ============================================================
%  LATENT HEAT — only the NH3 that condenses (2941.18 kmol/day)
%  2.941 kmol/day stays as vapor and does NOT contribute latent heat
%% ============================================================

F_NH3_cond_s = F_NH3_condensed / 86400;        % kmol/s
Q_latent     = -F_NH3_cond_s * delta_H_vap;    % kW (negative = heat removed)

%% ============================================================
%  TOTAL CONDENSER DUTY
%% ============================================================

Q_total_kW  = Q_sensible_total + Q_latent;
Q_total_MW  = Q_total_kW / 1000;
Q_total_MJd = Q_total_kW * 86.4;

%% ============================================================
%  DISPLAY RESULTS
%% ============================================================

fprintf('=======================================================\n');
fprintf('              CONDENSER DUTY RESULTS                  \n');
fprintf('=======================================================\n');
fprintf('Condenser pressure         : %7.1f bar  (< NH3 Pc=113.5 bar)\n', P_cond);
fprintf('Inlet temperature  (T_in)  : %7.2f K  (%6.2f deg C)\n', T_in,  T_in-273.15);
fprintf('-------------------------------------------------------\n');
fprintf('NH3 vapor split specification:\n');
fprintf('  NH3 condensed (liquid)   : %7.2f kmol/day  -> product\n', F_NH3_condensed);
fprintf('  NH3 vapor (recycle gas)  : %7.3f kmol/day  -> recycle\n', F_NH3_vapor);
fprintf('  NH3 mole frac in gas     : %.4e\n', y_NH3_gas);
fprintf('  NH3 partial pressure     : %.4e bar  (%.4f mmHg)\n', P_NH3_bar, P_NH3_mmHg);
fprintf('-------------------------------------------------------\n');
fprintf('Antoine back-calculation:\n');
fprintf('  T_out = B/(A - log10(P_mmHg)) - C\n');
fprintf('        = %.2f/(%.4f - log10(%.4f)) - %.2f\n', Ant_B, Ant_A, P_NH3_mmHg, Ant_C);
fprintf('Outlet temperature (T_out) : %7.2f K  (%6.2f deg C)\n', T_out, T_out_C);
fprintf('  Note: Deep refrigeration required at this temperature.\n');
fprintf('-------------------------------------------------------\n');
fprintf('ABCD coefficients (Cp/R = A+BT+CT^2+D/T^2):\n');
fprintf('  NH3: A=%7.4f  B=%.4e  C=%.4e  D=%.4e\n', ABCD_NH3);
fprintf('  N2 : A=%7.4f  B=%.4e  C=%.4e  D=%.4e\n', ABCD_N2);
fprintf('  H2 : A=%7.4f  B=%.4e  C=%.4e  D=%.4e\n', ABCD_H2);
fprintf('-------------------------------------------------------\n');
fprintf('dH sensible - NH3          : %10.2f kJ/kmol\n', dH_NH3);
fprintf('dH sensible - N2           : %10.2f kJ/kmol\n', dH_N2);
fprintf('dH sensible - H2           : %10.2f kJ/kmol\n', dH_H2);
fprintf('-------------------------------------------------------\n');
fprintf('Sensible duty - NH3        : %10.2f kW\n', Q_sensible_NH3);
fprintf('Sensible duty - N2         : %10.2f kW\n', Q_sensible_N2);
fprintf('Sensible duty - H2         : %10.2f kW\n', Q_sensible_H2);
fprintf('Total sensible duty        : %10.2f kW\n', Q_sensible_total);
fprintf('-------------------------------------------------------\n');
fprintf('delta H vap NH3 at 80 bar  : %10.2f kJ/kmol\n', delta_H_vap);
fprintf('NH3 condensed              : %10.2f kmol/day\n', F_NH3_condensed);
fprintf('Latent heat duty           : %10.2f kW\n', Q_latent);
fprintf('-------------------------------------------------------\n');
fprintf('TOTAL CONDENSER DUTY       : %10.2f kW\n',  Q_total_kW);
fprintf('  (negative = heat removed): %10.4f MW\n',  Q_total_MW);
fprintf('                           : %10.2f MJ/day\n', Q_total_MJd);
fprintf('=======================================================\n');
fprintf('\nFeed to separator: %.2f deg C, %.0f bar\n', T_out_C, P_cond);