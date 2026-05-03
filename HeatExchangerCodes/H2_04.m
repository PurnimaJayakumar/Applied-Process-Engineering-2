%% ---------- ENERGY BALANCE: STEAM METHANE REFORMING (SMR) ----------
clear; clc;

%% -------------------- Cp/R COEFFICIENTS (A,B,C,D) --------------------
% Cp/R = A + B*T + C*T^2 + D/T^2
%
coeffs_CH4 = {1.702, 9.081E-3, -2.164E-6, 0};
coeffs_H2O = {3.470, 1.450E-3, 0, 0.121E+5};
coeffs_CO  = {3.376, 0.557E-3, 0, -0.031E+5};
coeffs_H2  = {3.249, 0.422E-3, 0, 0.083E+5};

R_val = 8.314;   % J/mol-K  (8.314 J/mol-K)

%% -------------------- INPUT MOLAR FLOW RATE (kmol/day) --------------------
n_CH4_in = 1102.94115;
n_H2O_in = 2205.8823;

%% Reactor conditions
xi = 1102.94115;         % extent of reaction
T_in  = 550 + 273.15;    % K
T_rxn = 900 + 273.15;    % K
T_ref = 298.15;          % reference

%% -------------------- OUTLET MOLAR FLOW RATE --------------------
% CH4 + H2O → CO + 3H2
n_CH4_out = n_CH4_in - xi;
n_H2O_out = n_H2O_in - xi;
n_CO_out  = xi;
n_H2_out  = 3*xi;

%% -------------------- SENSIBLE HEAT OF FEED (T_in → T_rxn) --------------------
Q_CH4_in = n_CH4_in * H_change(T_in, T_rxn, coeffs_CH4{:}, R_val);
Q_H2O_in = n_H2O_in * H_change(T_in, T_rxn, coeffs_H2O{:}, R_val);

Q_sensible_in = Q_CH4_in + Q_H2O_in;

%% -------------------- REACTION HEAT @ T_rxn --------------------
% Base ΔH at 298 K
DeltaH_SMR_298 = +206.1;     % kJ/mol

% Temperature-dependent enthalpy increments from 298 → T_rxn
dH_CH4 = H_change(T_ref, T_rxn, coeffs_CH4{:}, R_val);
dH_H2O = H_change(T_ref, T_rxn, coeffs_H2O{:}, R_val);
dH_CO  = H_change(T_ref, T_rxn, coeffs_CO{:},  R_val);
dH_H2  = H_change(T_ref, T_rxn, coeffs_H2{:},  R_val);

% ΔH(T) = ΔH(298) + Σν_i Δh_i(T)
DeltaH_T = DeltaH_SMR_298 + (1*dH_CO + 3*dH_H2) - (1*dH_CH4 + 1*dH_H2O);

% Total reaction heat at reactor temperature
Q_rxn = DeltaH_T * xi;

%% -------------------- TOTAL ENERGY DUTY --------------------
Q_total = Q_sensible_in + Q_rxn;

%% -------------------- DISPLAY --------------------
disp("----- SMR ENERGY BALANCE RESULTS -----");
fprintf('Feed sensible heat (T_in → T_rxn):  %E kJ/day\n', Q_sensible_in);
fprintf('Reaction heat ΔH_rxn(1173 K):       %E kJ/day\n', Q_rxn);
fprintf('ΔH_rxn at 1173 K (per mol):         %f kJ/mol\n', DeltaH_T);
fprintf('-------------------------------------------------\n');
fprintf('TOTAL DUTY REQUIRED:                %E kJ/day\n', Q_total);


%% ----------- FUNCTION: ENTHALPY CHANGE BY Cp INTEGRAL -------------
function delta_H = H_change(T1, T2, A, B, C, D, R_val)
    % Cp/R = A + B*T + C*T^2 + D/T^2
    % H(T)/R = ∫(Cp/R)dT = A*T + B*T^2/2 + C*T^3/3 - D/T

    H_R_T2 = A*T2 + B*T2^2/2 + C*T2^3/3 - D/T2;
    H_R_T1 = A*T1 + B*T1^2/2 + C*T1^3/3 - D/T1;

    delta_H_over_R = H_R_T2 - H_R_T1;

    delta_H = delta_H_over_R * R_val;  % kJ/mol
end
