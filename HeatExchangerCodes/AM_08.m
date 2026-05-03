clear; clc;

%% ---------- CONSTANTS ----------
R = 8.314;  % kJ/kmol-K  (for Cp integrals)



%% ---------- Cp/R POLYNOMIAL COEFFICIENTS ----------
% Cp/R = A + B T + C T^2 + D/T^2
coeffs_N2   = {3.280, 0.593E-3, 0, 0.040E+5};
coeffs_H2   = {3.249, 0.422E-3, 0, 0.083E+5};
coeffs_NH3g = {3.578, 3.020E-3, 0, -0.186E+5};


%% ---------- USER INPUTS ----------
% kmol/day

nN2_fresh = 5882.3528;
nH2_fresh = 17647.0584;
nNH3_fresh = 2.941;
% Temperatures (K)
T_in = 245.67;
T_out = 450+273;

Q_sens_N2_fresh = nN2_fresh * H_change(T_in, T_out, coeffs_N2{:}, R);
Q_sens_H2_fresh = nH2_fresh * H_change(T_in, T_out, coeffs_H2{:}, R);
Q_sens_NH3_fresh = nNH3_fresh * H_change(T_in, T_out, coeffs_NH3g{:}, R);


fprintf('Q_sens N2 fresh      : %E kJ/day\n', Q_sens_N2_fresh);
fprintf('Q_sens H2 fresh      : %E kJ/day\n', Q_sens_H2_fresh);
fprintf('Q_sens NH3 fresh      : %E kJ/day\n', Q_sens_NH3_fresh);
fprintf('Q_sens total         : %E kJ/day\n', Q_sens_H2_fresh+Q_sens_N2_fresh+Q_sens_NH3_fresh);
function delta_H = H_change(T1, T2, A, B, C, D, R_val)
    H_over_R_T2 = A*T2 + B*T2^2/2 + C*T2^3/3 - D/T2;
    H_over_R_T1 = A*T1 + B*T1^2/2 + C*T1^3/3 - D/T1;
    delta_H = (H_over_R_T2 - H_over_R_T1) * R_val;
end
