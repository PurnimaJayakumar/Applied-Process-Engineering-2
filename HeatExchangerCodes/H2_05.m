%% WGS_FULL_ENERGY_BALANCE.M
% Cp/R = A + B*T + C*T^2 + D/T^2
% Cp   = R * (A + B*T + C*T^2 + D/T^2)
% Units consistent: R in J/mol/K, flows in kmol/day, output in kJ/day.

clear; clc;

%% -------------------- INPUT DATA --------------------
n_H2_in  = 3308.82345;     % kmol/day
n_H2O_in = 1102.94115;
n_CO_in  = 1102.94115;
n_CO2_in = 0;
n_CH4_in = 0;

xi_kmol_per_day = 1102.94115;   % kmol/day reacted

T_in  = 350 + 273.15;  % K
T_rxn = 700;           % K  (effective)
T_out = 260 + 273.15;  % K
T0    = 298.15;        % K reference

DeltaH_rxn_298 = -41.1;  % kJ/mol  (CO + H2O -> CO2 + H2)

%% -------------------- COEFFICIENTS FOR Cp/R --------------------
coeff.H2  = [3.249, 0.422e-3, 0,        0.083e5];
coeff.H2O = [3.470, 1.450e-3, 0,        0.121e5];
coeff.CO  = [3.376, 0.557e-3, 0,       -0.031e5];
coeff.CO2 = [5.457, 1.045e-3, 0,       -1.157e5];
coeff.CH4 = [1.702, 9.081e-3, -2.164e-6, 0];

R = 8.314;   % J/mol/K

%% -------------------- FUNCTIONS --------------------
% Cp(T) [J/mol/K]
Cp = @(c, T) R * ( c(1) + c(2)*T + c(3)*T.^2 + c(4)./T.^2 );

% Integral Cp dT (J/mol)
delta_h = @(c, T1, T2) R * (c(1)*(T2-T1) + 0.5*c(2)*(T2^2 - T1^2) +(1/3)*c(3)*(T2^3 - T1^3) -c(4)*(1/T2 - 1/T1) );

%% -------------------- OUTLET FLOWS --------------------
n_CO_out  = n_CO_in  - xi_kmol_per_day;
n_H2O_out = n_H2O_in - xi_kmol_per_day;
n_CO2_out = n_CO2_in + xi_kmol_per_day;
n_H2_out  = n_H2_in  + xi_kmol_per_day;
n_CH4_out = n_CH4_in;

%% -------------------- FEED HEATING: T_in → T_rxn --------------------
Q_sensible_in = ...
    (n_H2_in  *1000)*delta_h(coeff.H2,  T_in, T_rxn) + ...
    (n_H2O_in *1000)*delta_h(coeff.H2O, T_in, T_rxn) + ...
    (n_CO_in  *1000)*delta_h(coeff.CO,  T_in, T_rxn) + ...
    (n_CO2_in *1000)*delta_h(coeff.CO2, T_in, T_rxn) + ...
    (n_CH4_in *1000)*delta_h(coeff.CH4, T_in, T_rxn);

Q_sensible_in = Q_sensible_in / 1000;   % → kJ/day

%% -------------------- REACTION HEAT AT T_rxn --------------------
dH_H2_T  = delta_h(coeff.H2,  T0, T_rxn)/1000;   % kJ/mol
dH_H2O_T = delta_h(coeff.H2O, T0, T_rxn)/1000;
dH_CO_T  = delta_h(coeff.CO,  T0, T_rxn)/1000;
dH_CO2_T = delta_h(coeff.CO2, T0, T_rxn)/1000;

DeltaH_rxn_T = DeltaH_rxn_298 + (dH_CO2_T + dH_H2_T) - (dH_CO_T + dH_H2O_T);

Q_rxn = DeltaH_rxn_T * (xi_kmol_per_day*1000);   % kJ/day

%% -------------------- COOLING PRODUCTS: T_rxn → T_out --------------------
Q_sensible_out = ...
    (n_H2_out  *1000)*delta_h(coeff.H2,  T_out, T_rxn) + ...
    (n_H2O_out *1000)*delta_h(coeff.H2O, T_out, T_rxn) + ...
    (n_CO_out  *1000)*delta_h(coeff.CO,  T_out, T_rxn) + ...
    (n_CO2_out *1000)*delta_h(coeff.CO2, T_out, T_rxn) + ...
    (n_CH4_out *1000)*delta_h(coeff.CH4, T_out, T_rxn);

Q_sensible_out = Q_sensible_out / 1000;   % kJ/day

%% -------------------- TOTAL --------------------
Q_total = Q_sensible_in + Q_rxn - Q_sensible_out;

%% -------------------- DISPLAY --------------------
fprintf('----- ENERGY BALANCE (WGS, Correct Units) -----\n');
fprintf('Feed sensible:           %+E kJ/day\n', Q_sensible_in);
fprintf('Reaction heat (@700 K):  %+E kJ/day   (ΔH = %.3f kJ/mol)\n', Q_rxn, DeltaH_rxn_T);
fprintf('Cooling products:        %+E kJ/day\n', Q_sensible_out);
fprintf('-----------------------------------------------\n');
fprintf('NET HEAT RELEASED:       %+E kJ/day\n', Q_total);
fprintf('-----------------------------------------------\n');
