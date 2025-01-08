function [V, N] = solveBenzenoODE(N0, V_span, T)
    % Define the system of ODEs
    function dNdV = benzenoODE(V, N)
        % Constants
        _k1 = 7.0e5; % L/mol.h
        K1 = 0.31;
        _k2 = 4.0e5; % L/mol.h
        K2 = 0.48;
        P = 101325; % Pa
        R = 8314;   % J/(mol.K)
        
        % Total moles
        N_total = sum(N);
        
        % Concentrations (L/mol)
        C_c6h6 = (P / (R * T)) * (N(1) / N_total);
        C_h = (P / (R * T)) * (N(2) / N_total);
        C_c12h10 = (P / (R * T)) * (N(3) / N_total);
        C_c18h14 = (P / (R * T)) * (N(4) / N_total);
        
        % Reaction rates
        r1 = _k1 * (C_c6h6^2 - (C_c12h10 * C_h) / K1);
        r2 = _k2 * (C_c6h6 * C_c12h10 - (C_c18h14 * C_h) / K2);
        
        % ODEs
        dNdV = zeros(4, 1);
        dNdV(1) = -(2 * r1) - r2;  % dNc6h6/dV
        dNdV(2) = r1 + r2;         % dNh/dV
        dNdV(3) = r1 - r2;         % dNc12h10/dV
        dNdV(4) = r2;              % dNc18h14/dV
    end

    % Solve the ODE
    [V, N] = ode45(@benzenoODE, V_span, N0);
end