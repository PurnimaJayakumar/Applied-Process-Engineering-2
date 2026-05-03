clear; clc;

%% ---------- CONSTANTS ----------
R = 8.314;  % kJ/kmol-K  (for Cp integrals)

% Reaction enthalpy for NH3 formation
dH_rxn_per_mol_NH3 = -46.2;  % kJ/mol NH3

%% ---------- Cp/R POLYNOMIAL COEFFICIENTS ----------
% Cp/R = A + B T + C T^2 + D/T^2
coeffs_N2   = {3.280, 0.593E-3, 0, 0.040E+5};
coeffs_H2   = {3.249, 0.422E-3, 0, 0.083E+5};
coeffs_NH3g = {3.578, 3.020E-3, 0, -0.186E+5};


%% ---------- USER INPUTS ----------
% kmol/day
T_reactor = 773;

nN2_fresh = 1470.5882;
nH2_fresh = 4411.7646;
% Temperatures (K)
T_in = 703.44;
T_out = T_reactor;

Q_sens_N2_fresh = nN2_fresh * H_change(T_in, T_out, coeffs_N2{:}, R);
Q_sens_H2_fresh = nH2_fresh * H_change(T_in, T_out, coeffs_H2{:}, R);

fprintf('-- Inlet sensible heats to mixing unit (T = %.1f K) --\n', T_reactor);
fprintf('Q_sens N2 fresh      : %E kJ/day\n', Q_sens_N2_fresh);
fprintf('Q_sens H2 fresh      : %E kJ/day\n', Q_sens_H2_fresh);
fprintf('Q_sens total         : %E kJ/day\n', Q_sens_H2_fresh+Q_sens_N2_fresh);
function delta_H = H_change(T1, T2, A, B, C, D, R_val)
    H_over_R_T2 = A*T2 + B*T2^2/2 + C*T2^3/3 - D/T2;
    H_over_R_T1 = A*T1 + B*T1^2/2 + C*T1^3/3 - D/T1;
    delta_H = (H_over_R_T2 - H_over_R_T1) * R_val;
end
