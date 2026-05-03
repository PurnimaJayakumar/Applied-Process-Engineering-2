%% =========================================================
%  WATER CONDENSATION FROM COMPRESSED AIR
%  Cooling compressed air from 494.95 K to 303.15 K (30 degC)
%  Operating pressure: 6 bar
%  Condensate: water vapour -> liquid water
%
%  Steps:
%  1. Check how much water remains as vapour at 30 degC, 6 bar
%     using Antoine equation (saturation pressure)
%  2. Calculate condensate amount
%  3. Calculate total heat duty Q (sensible cooling + latent heat)
%  4. Report exit stream composition on dry basis
%% =========================================================

clc; clear; close all;
R = 8.314;  % J/(mol.K)

%% -------------------------------------------------------
%  1. MOLAR FLOW RATES  [kmol/day -> mol/day]
% -------------------------------------------------------
n_N2  = 1472.06026e3;      % mol/day
n_O2  = 390.9959761e3;     % mol/day
n_CO2 = 0.746549694e3;     % mol/day
n_H2O = 35.25373555e3;     % mol/day

n_total_in = n_N2 + n_O2 + n_CO2 + n_H2O;  % mol/day

y_N2  = n_N2  / n_total_in;
y_O2  = n_O2  / n_total_in;
y_CO2 = n_CO2 / n_total_in;
y_H2O = n_H2O / n_total_in;

T_in  = 494.95;   % K  inlet to condenser (exit of compressor)
T_out = 303.15;   % K  = 30 degC  target outlet temperature
P_tot = 6;        % bar  operating pressure

fprintf('============================================================\n');
fprintf('  INLET STREAM TO CONDENSER\n');
fprintf('============================================================\n');
fprintf('T_in  = %.2f K  (%.2f degC)\n', T_in,  T_in-273.15);
fprintf('T_out = %.2f K  (%.2f degC)\n', T_out, T_out-273.15);
fprintf('P     = %.1f bar\n\n', P_tot);
fprintf('n_N2  = %.4f kmol/day\n', n_N2/1e3);
fprintf('n_O2  = %.4f kmol/day\n', n_O2/1e3);
fprintf('n_CO2 = %.4f kmol/day\n', n_CO2/1e3);
fprintf('n_H2O = %.4f kmol/day\n', n_H2O/1e3);
fprintf('n_total_in = %.4f kmol/day\n\n', n_total_in/1e3);

%% -------------------------------------------------------
%  2. SATURATION PRESSURE OF WATER AT 30 degC
%    
%     ln10(P_sat/bar) = A - B/(C + T_Cels)
%    
%    
%     

% -------------------------------------------------------
A_ant = 16.3872;
B_ant = 3885.70;
C_ant = 230.170;

P_sat_water = exp(A_ant - B_ant / (C_ant + (T_out-273.15) )   )*1e-2;  % bar

fprintf('============================================================\n');
fprintf('  PHASE EQUILIBRIUM AT 30 degC, 6 bar\n');
fprintf('============================================================\n');
fprintf('Saturation pressure of water at 30 degC = %.4f bar\n', P_sat_water);

% Partial pressure of water in inlet stream
P_H2O_in = y_H2O * P_tot;  % bar
fprintf('Partial pressure of H2O in inlet stream = %.4f bar\n\n', P_H2O_in);

%% -------------------------------------------------------
%  3. CONDENSATION CALCULATION
%     At equilibrium, vapour phase partial pressure of water
%     cannot exceed P_sat at T_out.
%
%     Max moles of H2O that can stay as vapour:
%     y_H2O_out = P_sat / P_tot  (Raoult / Dalton)
%
%     Let n_H2O_vap = moles remaining as vapour
%     Let n_dry = total moles of non-condensable gases (N2+O2+CO2)
%
%     y_H2O_out = n_H2O_vap / (n_dry + n_H2O_vap) = P_sat/P_tot
%     Solving: n_H2O_vap = n_dry * (P_sat/P_tot) / (1 - P_sat/P_tot)
% -------------------------------------------------------
n_dry = n_N2 + n_O2 + n_CO2;  % mol/day  (non-condensables)

y_H2O_sat = P_sat_water / P_tot;  % max mole fraction of water vapour at outlet

% Moles of water remaining as vapour
n_H2O_vap = n_dry * y_H2O_sat / (1 - y_H2O_sat);  % mol/day

% Moles condensed
n_H2O_cond = n_H2O - n_H2O_vap;  % mol/day

% Check: if inlet water is already below saturation, nothing condenses
if n_H2O_cond < 0
    fprintf('Water partial pressure (%.4f bar) < P_sat (%.4f bar)\n', P_H2O_in, P_sat_water);
    fprintf('No condensation occurs at these conditions.\n\n');
    n_H2O_cond = 0;
    n_H2O_vap  = n_H2O;
else
    fprintf('Condensation occurs: P_H2O_in (%.4f bar) > P_sat (%.4f bar)\n\n', ...
             P_H2O_in, P_sat_water);
end

fprintf('============================================================\n');
fprintf('  CONDENSATION RESULTS\n');
fprintf('============================================================\n');
fprintf('H2O in inlet stream      = %.4f kmol/day\n', n_H2O/1e3);
fprintf('H2O remaining as vapour  = %.4f kmol/day\n', n_H2O_vap/1e3);
fprintf('H2O condensed (liquid)   = %.4f kmol/day\n', n_H2O_cond/1e3);
fprintf('H2O condensed (kg/day)   = %.4f kg/day\n\n', n_H2O_cond/1e3 * 18);

%% -------------------------------------------------------
%  4. SHOMATE COEFFICIENTS FOR Cp
%     Cp (J/mol.K) = A + B*T + C*T^2 + D/T^2
% -------------------------------------------------------
coeff_N2  = R*[3.280,  0.593e-3,  0,  0.040e+5];
coeff_O2  = R*[3.639,  0.506e-3,  0, -0.227e+5];
coeff_CO2 = R*[5.457,  1.045e-3,  0, -1.157e+5];
coeff_H2O = R*[3.470,  1.450e-3,  0,  0.121e+5];  % vapour phase Cp

Cp_N2  = @(T) coeff_N2(1)  + coeff_N2(2).*T  + coeff_N2(3).*T.^2  + coeff_N2(4)./T.^2;
Cp_O2  = @(T) coeff_O2(1)  + coeff_O2(2).*T  + coeff_O2(3).*T.^2  + coeff_O2(4)./T.^2;
Cp_CO2 = @(T) coeff_CO2(1) + coeff_CO2(2).*T + coeff_CO2(3).*T.^2 + coeff_CO2(4)./T.^2;
Cp_H2O = @(T) coeff_H2O(1) + coeff_H2O(2).*T + coeff_H2O(3).*T.^2 + coeff_H2O(4)./T.^2;

%% -------------------------------------------------------
%  5. HEAT DUTY CALCULATION
%
%  Three contributions:
%  (a) Sensible cooling of ALL gas from T_in to T_out
%      Q_sensible = sum[ n_i * integral(T_in->T_out) Cp_i dT ]
%      (for H2O, only the vapour portion n_H2O_vap cools as gas
%       the condensing portion n_H2O_cond cools to T_out then condenses)
%
%  (b) Latent heat of condensation of n_H2O_cond
%      Q_latent = n_H2O_cond * delta_Hvap(at T_out)
%
%  delta_Hvap of water at 30 degC = 2430.7 J/g = 2430.7 * 18 J/mol
%  (from steam tables at 30 degC)
% -------------------------------------------------------
dH_vap_water = 2430.7 * 18;  % J/mol  latent heat at 30 degC

% Sensible cooling: all components cool from T_in to T_out
Q_N2  = n_N2       * integral(@(T) Cp_N2(T),  T_in, T_out);  % J/day  (negative = heat removed)
Q_O2  = n_O2       * integral(@(T) Cp_O2(T),  T_in, T_out);
Q_CO2 = n_CO2      * integral(@(T) Cp_CO2(T), T_in, T_out);
Q_H2O_sens = n_H2O * integral(@(T) Cp_H2O(T), T_in, T_out);  % all H2O cools first

% Latent heat released during condensation (negative = heat removed)
Q_latent = -n_H2O_cond * dH_vap_water;  % J/day

% Total heat duty (negative = removed from stream, i.e. cooling)
Q_sensible_total = Q_N2 + Q_O2 + Q_CO2 + Q_H2O_sens;  % J/day
Q_total = Q_sensible_total + Q_latent;                  % J/day

fprintf('============================================================\n');
fprintf('  HEAT DUTY\n');
fprintf('============================================================\n');
fprintf('Sensible cooling duty:\n');
fprintf('  Q_N2       = %.4e  kJ/day\n', Q_N2/1e3);
fprintf('  Q_O2       = %.4e  kJ/day\n', Q_O2/1e3);
fprintf('  Q_CO2      = %.4e  kJ/day\n', Q_CO2/1e3);
fprintf('  Q_H2O sens = %.4e  kJ/day\n', Q_H2O_sens/1e3);
fprintf('  Q_sensible = %.4e  kJ/day\n\n', Q_sensible_total/1e3);
fprintf('Latent heat of condensation:\n');
fprintf('  n_H2O condensed    = %.4f kmol/day\n', n_H2O_cond/1e3);
fprintf('  dH_vap at 30 degC  = %.2f J/mol\n', dH_vap_water);
fprintf('  Q_latent           = %.4e  kJ/day\n\n', Q_latent/1e3);
fprintf('TOTAL heat duty Q  = %.4e  kJ/day\n',   Q_total/1e3);
fprintf('                   = %.4e  MJ/day\n',   Q_total/1e6);
fprintf('(negative = heat removed from stream)\n\n');

% As power
Q_power_kW = abs(Q_total) / 86400 / 1e3;  % kW
fprintf('Cooling power required = %.4f kW = %.6f MW\n\n', Q_power_kW, Q_power_kW/1e3);

%% -------------------------------------------------------
%  6. EXIT STREAM COMPOSITION
% -------------------------------------------------------
n_total_out = n_N2 + n_O2 + n_CO2 + n_H2O_vap;  % gas phase only

y_N2_out  = n_N2       / n_total_out;
y_O2_out  = n_O2       / n_total_out;
y_CO2_out = n_CO2      / n_total_out;
y_H2O_out = n_H2O_vap  / n_total_out;

fprintf('============================================================\n');
fprintf('  EXIT STREAM COMPOSITION (gas phase)\n');
fprintf('============================================================\n');
fprintf('%-10s %15s %15s %15s\n', 'Component', 'kmol/day', 'Mole frac in', 'Mole frac out');
fprintf('%s\n', repmat('-',1,58));
fprintf('%-10s %15.4f %15.4f %15.4f\n', 'N2',  n_N2/1e3,       y_N2,  y_N2_out);
fprintf('%-10s %15.4f %15.4f %15.4f\n', 'O2',  n_O2/1e3,       y_O2,  y_O2_out);
fprintf('%-10s %15.4f %15.4f %15.4f\n', 'CO2', n_CO2/1e3,      y_CO2, y_CO2_out);
fprintf('%-10s %15.4f %15.4f %15.4f\n', 'H2O', n_H2O_vap/1e3,  y_H2O, y_H2O_out);
fprintf('%s\n', repmat('-',1,58));
fprintf('%-10s %15.4f\n', 'Condensate', n_H2O_cond/1e3);
fprintf('  (liquid water removed = %.4f kg/day)\n\n', n_H2O_cond/1e3*18);

%% -------------------------------------------------------
%  7. SUMMARY TABLE
% -------------------------------------------------------
fprintf('============================================================\n');
fprintf('  SUMMARY\n');
fprintf('============================================================\n');
fprintf('%-40s %15s\n', 'Quantity', 'Value');
fprintf('%s\n', repmat('-',1,57));
fprintf('%-40s %15.2f  K\n',      'T inlet',                  T_in);
fprintf('%-40s %15.2f  K\n',      'T outlet',                 T_out);
fprintf('%-40s %15.2f  bar\n',    'Pressure',                 P_tot);
fprintf('%-40s %15.4f  bar\n',    'P_sat water at 30 degC',   P_sat_water);
fprintf('%-40s %15.4f  kmol/day','H2O condensed',             n_H2O_cond/1e3);
fprintf('\n');
fprintf('%-40s %15.4f  kg/day\n', 'Condensate mass',          n_H2O_cond/1e3*18);
fprintf('%-40s %15.4e  kJ/day\n', 'Total heat duty',          Q_total/1e3);
fprintf('%-40s %15.4f  kW\n',     'Cooling power',            Q_power_kW);
fprintf('%-40s %15.6f  MW\n',     'Cooling power',            Q_power_kW/1e3);