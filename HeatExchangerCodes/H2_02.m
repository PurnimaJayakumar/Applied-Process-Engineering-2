clear; clc;

%% ---------- CONSTANTS ----------
R = 8.314;   % kJ/kmol-K  (for Cp integrals)


%% ---------- Cp/R POLYNOMIAL COEFFICIENTS ----------
% Cp/R = A + B T + C T^2 + D/T^2
coeffs_CH4   = {1.702, 9.081E-3, -2.164E-6, 0};
coeffs_H2O  = {3.470,     1.450E-3,    0,        0.121E+5};

%% ---------- USER INPUTS ----------
% kmol/day

nCH4 = 1102.94115; 
nH2O = 2205.8823;
% Temperatures (K)
T_in = 673; % 400 deg Cels
T_out = 823; % 550 deg Cels

Q_sens_CH4 = nCH4 * H_change(T_in, T_out, coeffs_CH4{:}, R);
Q_sens_H2O = nH2O * H_change(T_in, T_out, coeffs_H2O{:}, R);

fprintf('-- Total (T = %.1f K) --\n', T_in);
fprintf('Q_sens total      : %E kJ/day\n', Q_sens_CH4+Q_sens_H2O);


function delta_H = H_change(T1, T2, A, B, C, D, R_val)
    H_over_R_T2 = A*T2 + B*T2^2/2 + C*T2^3/3 - D/T2;
    H_over_R_T1 = A*T1 + B*T1^2/2 + C*T1^3/3 - D/T1;
    delta_H = (H_over_R_T2 - H_over_R_T1) * R_val;
end
