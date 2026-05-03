clear; clc;

%% ---------- CONSTANTS ----------
R = 8.314;   % kJ/kmol-K  (for Cp integrals)


%% ---------- Cp/R POLYNOMIAL COEFFICIENTS ----------
% Cp/R = A + B T + C T^2 + D/T^2
coeffs_H2  = {3.249, 0.422E-3, 0, 0.083E+5};

%% ---------- USER INPUTS ----------
% kmol/day

nH2 = 4411.7646;
% Temperatures (K)
T_in = 260+273; % 400 deg Cels
T_out = 27+273; % 550 deg Cels

Q_sens_H2 = nH2 * H_change(T_in, T_out, coeffs_H2{:}, R);

fprintf('-- Total (T = %.1f K) --\n', T_in);
fprintf('Q_sens total      : %E kJ/day\n', Q_sens_H2);


function delta_H = H_change(T1, T2, A, B, C, D, R_val)
    H_over_R_T2 = A*T2 + B*T2^2/2 + C*T2^3/3 - D/T2;
    H_over_R_T1 = A*T1 + B*T1^2/2 + C*T1^3/3 - D/T1;
    delta_H = (H_over_R_T2 - H_over_R_T1) * R_val;
end
