clc; clear all;

R_val = 8.314; % J/mol-K 

%% ---------------- Cp COEFFICIENTS (A, B, C, D) ----------------
 
coeffs_H2   = {3.249,     0.422E-3,    0,        0.083E+5};      
coeffs_CO   = {3.376,     0.557E-3,    0,       -0.031E+5};
coeffs_H2O  = {3.470,     1.450E-3,    0,        0.121E+5};
T1= 900+273;%K
T2=350+273; %K

nH2  = 3308.82345;
nH2O = 1102.94115;
nCO = 1102.94115;

H_change = @(T1,T2,A,B,C,D,R) ...
    R*( A*(T2-T1) ...
       + B*(0.5*(T2^2 - T1^2)) ...
       + C*(1/3*(T2^3 - T1^3)) ...
       + D*(-1*(1./T2 - 1./T1)) );
Q_sens_1 = nCO * H_change(T1,T2,coeffs_CO{:},R_val);
Q_sens_2 = nH2O * H_change(T1,T2,coeffs_H2O{:},R_val);
Q_sens_3 = nH2 * H_change(T1,T2,coeffs_H2{:},R_val);

Q_total = Q_sens_1 + Q_sens_2 + Q_sens_3;

fprintf('Heat Exchanger (Hydrogen: After SMR Before WGS): %E kJ/day \n', Q_total);