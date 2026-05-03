% Cryo Energy balance
%% =========================================================
%  CRYO UNIT HEAT BALANCE
%  Feed: Air (N2 + O2) at -170 deg C, 6 bar
%  Outlet: N2 gas + O2 liquid at -183 deg C
%  Fill in your heat values below
%% =========================================================

clc; clear; close all;

%% -------------------------------------------------------
%  1. MOLAR FLOW RATES  [kmol/day]
% -------------------------------------------------------
n_N2 = 1472.06026;    % kmol/day
n_O2 = 390.9959761;   % kmol/day

%% -------------------------------------------------------
%  2. TEMPERATURES
% -------------------------------------------------------
T_in  = 103.15;   T1 = T_in;
T_out = 90.15;    T2 = T_out;

%% -------------------------------------------------------
%  3. FILL IN YOUR HEAT VALUES  [kJ/kmol]
%     All values should be NEGATIVE (heat removed)
% -------------------------------------------------------
Cp_N2 = @(T) 28.98641 + 1.853978*(T/1000) + (-9.647459)*(T/1000).^2 + 16.63537*(T/1000).^3 + 0.000117./(T/1000).^2;
Cp_O2 = @(T) 31.32234	 + -20.23531*(T/1000) + 57.86644*(T/1000).^2 + (-36.50624)*(T/1000).^3 + (-0.007374)./(T/1000).^2;
dH_sens_N2 = integral(Cp_N2, T1, T2);  % kJ/kmol
dH_sens_O2 = integral(Cp_O2, T1, T2);  % kJ/kmol
      
%LATENT HEAT OF VAPORIZATION PER KG O2: 213000 J/kg 
% = 213000 Jkg * 32E-3 kg/mol = 6816 J/mol 
% Src:https://www.engineeringtoolbox.com/oxygen-d_1422.html
% --- O2 latent heat of liquefaction at -183 degC ---

dH_lat_O2  = -6816;        % kJ/kmol 

%% -------------------------------------------------------
%  4. TOTAL HEAT DUTIES  [kJ/day]
% -------------------------------------------------------
Q_sens_N2  = dH_sens_N2 * n_N2;    % kJ/day
Q_sens_O2  = dH_sens_O2 * n_O2;    % kJ/day
Q_lat_O2   = dH_lat_O2  * n_O2;    % kJ/day

Q_sens_total = Q_sens_N2 + Q_sens_O2;          % kJ/day
Q_lat_total  = Q_lat_O2;                        % kJ/day
Q_total      = Q_sens_total + Q_lat_total;      % kJ/day

%% -------------------------------------------------------
%  5. POWER
% -------------------------------------------------------
Q_power_kW = abs(Q_total) / 86400;   % kJ/day ÷ s/day = kW  ✓
Q_power_MW = Q_power_kW / 1000;      % MW

%% -------------------------------------------------------
%  6. OUTLET STREAM CONDITIONS
% -------------------------------------------------------
P_N2_out = 6;        % bar  (gaseous N2, stays at feed pressure)
P_O2_out = 1.013;    % bar  (liquid O2, flashes to atmospheric)
T_N2_out = 90.15;    % K = -183 degC
T_O2_out = 90.15;    % K = -183 degC

%% -------------------------------------------------------
%  7. RESULTS
% -------------------------------------------------------
fprintf('============================================================\n');
fprintf('  CRYO UNIT HEAT BALANCE\n');
fprintf('============================================================\n');
fprintf('Feed temperature  = %.2f K  (%.2f degC)\n', T_in,  T_in-273.15);
fprintf('Outlet temperature= %.2f K  (%.2f degC)\n', T_out, T_out-273.15);
fprintf('Delta T           = %.2f K\n\n', T_out - T_in);

fprintf('%-40s %15s\n', 'Quantity', 'Value');
fprintf('%s\n', repmat('-',1,57));
fprintf('%-40s %15.4e  kJ/day\n', 'Q sensible N2',    Q_sens_N2);
fprintf('%-40s %15.4e  kJ/day\n', 'Q sensible O2',    Q_sens_O2);
fprintf('%-40s %15.4e  kJ/day\n', 'Q sensible total', Q_sens_total);
fprintf('%-40s %15.4e  kJ/day\n', 'Q latent O2',      Q_lat_O2);
fprintf('%-40s %15.4e  kJ/day\n', 'Q TOTAL',          Q_total);
fprintf('%s\n', repmat('-',1,57));
fprintf('%-40s %15.4f  kW\n',     'Cooling power',    Q_power_kW);
fprintf('%-40s %15.6f  MW\n\n',   'Cooling power',    Q_power_MW);

fprintf('============================================================\n');
fprintf('  OUTLET STREAM SUMMARY\n');
fprintf('============================================================\n');
fprintf('%-12s %-10s %-12s %-12s %-15s\n', ...
        'Stream','Phase','Temp (K)','Temp (C)','Pressure (bar)');
fprintf('%s\n', repmat('-',1,63));
fprintf('%-12s %-10s %-12.2f %-12.2f %-15.3f\n', ...
        'N2','Gas',    T_N2_out, T_N2_out-273.15, P_N2_out);
fprintf('%-12s %-10s %-12.2f %-12.2f %-15.3f\n', ...
        'O2','Liquid', T_O2_out, T_O2_out-273.15, P_O2_out);