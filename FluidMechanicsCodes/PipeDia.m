function row = PipeDia(P_bar, T_C, n_kmol_day, y, MW, Tc_K, Pc_Pa, omega, L, Vmax, dP_allow_bar, eps, mu)
% =========================================================================
% PIPE DIAMETER CALCULATOR — Peng-Robinson EOS + Darcy-Weisbach
% =========================================================================
%
% INPUTS  (13 total)
% ------------------
%   P_bar         : Operating pressure                          (bar)
%   T_C           : Operating temperature                       (°C)
%   n_kmol_day    : Total molar flow rate                       (kmol/day)
%   y             : Mole fraction vector [y1, y2, ...]          (-)
%   MW            : Mixture molecular weight                    (kg/kmol)
%   Tc_K          : Critical temperatures, one per component    (K)
%   Pc_Pa         : Critical pressures,    one per component    (Pa)
%   omega         : Acentric factors,      one per component    (-)
%   L             : Pipe length                                 (m)
%   Vmax          : Maximum allowable velocity                  (m/s)
%   dP_allow_bar  : Maximum allowable pressure drop             (bar)
%   eps           : Absolute pipe roughness (commercial steel = 4.5e-5 m)
%   mu            : Dynamic viscosity                           (Pa.s)
%
% OUTPUTS (7 total)
% -----------------
%   D             : Pipe internal diameter                      (m)
%   Z             : PR compressibility factor                   (-)
%   V             : Fluid velocity                              (m/s)
%   dP_bar        : Pressure drop                               (bar)
%   Re            : Reynolds number                             (-)
%   rho           : Fluid density                               (kg/m3)
%   Q             : Volumetric flow rate                        (m3/s)
%
% COMPONENT PROPERTY REFERENCE
% -----------------------------
%   Species   Tc (K)    Pc (Pa)      omega
%   H2         33.2     12.80e5     -0.216
%   N2        126.2     33.90e5      0.037
%   O2        154.6     50.40e5      0.022
%   Ar        150.8     48.70e5      0.000
%   CO2       304.2     73.80e5      0.225
%   NH3       405.5    113.00e5      0.253
%
% EXAMPLE CALL
% ------------
%   Tc_K  = [126.2, 154.6];
%   Pc_Pa = [33.90, 50.40] * 1e5;
%   om    = [0.037, 0.022];
%   [D,Z,V,dP,Re] = PipeDia(6, 50, 1863, [0.79 0.21], 28.85, ...
%                            Tc_K, Pc_Pa, om, 10, 20, 0.3, 4.5e-5, 1.96e-5);
% =========================================================================

 
    nc = length(y);
   
    %% 1 — Unit conversions
    R        = 8.314;
    P        = P_bar * 1e5;               % Pa
    T        = T_C + 273.15;              % K
    n        = n_kmol_day * 1000 / 86400; % mol/s
    dP_allow = dP_allow_bar * 1e5;        % Pa

    %% 2 — Peng-Robinson EOS
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
    alpha = (1 + kappa .* (1 - sqrt(T ./ Tc_K))).^2;
    a_i   = 0.45724 * (R^2 * Tc_K.^2 ./ Pc_Pa);
    b_i   = 0.07780 * (R   * Tc_K    ./ Pc_Pa);
    aT_i  = a_i .* alpha;

    % van der Waals mixing rules
    a_mix = 0;
    for i = 1:nc
        for j = 1:nc
            a_mix = a_mix + y(i) * y(j) * sqrt(aT_i(i) * aT_i(j));
        end
    end
    b_mix = sum(y .* b_i);

    A_pr = a_mix * P / (R^2 * T^2);
    B_pr = b_mix * P / (R * T);

    coeff   = [1, -(1-B_pr), (A_pr - 3*B_pr^2 - 2*B_pr), -(A_pr*B_pr - B_pr^2 - B_pr^3)];
    Z_roots = roots(coeff);
    Z_real  = real(Z_roots(abs(imag(Z_roots)) < 1e-6 & real(Z_roots) > 0));

    if isempty(Z_real)
        error('PipeDia:PRfail', 'No valid Z root at P=%.1f bar, T=%.1f C.', P_bar, T_C);
    end
    Z = max(Z_real);

    %% 3 — Fluid properties
    rho = (P * MW) / (Z * R * T * 1000); % kg/m3
    Q   = n * Z * R * T / P;             % m3/s

    %% 4 — Minor losses
    K_total = 0.04 + 1.00;               % rounded entrance + pipe exit

    %% 5 — Iterative diameter search (1 mm steps)
    D_max     = 2.000;
    D_step    = 0.001;
    converged = false;
    D         = 0.010;

    while D <= D_max
        A_pipe = pi * D^2 / 4;
        V      = Q / A_pipe;
        Re     = rho * V * D / mu;

        % Colebrook-White friction factor (Darcy)
        if Re < 2300
            f = 64 / Re;
        else
            f = 0.02;
            for iter = 1:50
                f_new = 1 / (-2*log10(eps/(3.7*D) + 2.51/(Re*sqrt(f))))^2;
                if abs(f_new - f) < 1e-10, break; end
                f = f_new;
            end
            f = f_new;
        end

        % Darcy-Weisbach + minor losses
        dP_Pa = (f*(L/D) + K_total) * (rho*V^2/2);

        if (V <= Vmax) && (dP_Pa <= dP_allow)
            converged = true;
            break
        end
        D = D + D_step;
    end

    if ~converged
        error('PipeDia:noSolution', ...
              'No diameter up to %.1f m satisfies constraints.\nTry increasing Vmax (%.1f m/s) or dP_allow_bar (%.3f bar).', ...
              D_max, Vmax, dP_allow_bar);
    end

  
    %% 6 — Outputs
    dP_bar = dP_Pa / 1e5;

    % ---- RETURN AS ROW (for table use) ----
    row = [D, Z, V, dP_bar, Re, rho, Q];

end
