function [x, N, V, M] = StructuralCalcs(Vehicle_Data, Mach, q, F_T, m, I, x_cg, dmdx, h, G)

x = linspace(0, Vehicle_Data(7), 1000)';
dx = x(2);
x_rel_cg = x-x_cg;
m_sec = cumtrapz(x, dmdx, 1);

for i = 1:11
    stri = string(i-1);
    filename = strcat("CDDatafile", stri);

    [W_A, W_N, F_A, F_N] = AeroDistribution(Vehicle_Data, Mach, q, F_T, filename);

    a_cg_N = F_N/m;
    a_cg_A = F_A/m;

    M_A = W_N.*x_rel_cg;
    M_A_tot = trapz(x, M_A);
    alpha = M_A_tot/I;

    N(:, i) = -cumtrapz(x, W_A, 1) + a_cg_A*m_sec;
    V(:, i) = -cumtrapz(x,W_N, 1) + a_cg_N*m_sec + (alpha*cumtrapz(x, dmdx.*(x-x_cg), 1));
    M(:, i) = -cumtrapz(x,V(:, i));
end

% Flutter Calcs
C_r = Vehicle_Data(9);
C_t = Vehicle_Data(22);
%b = Vehicle_Data(25);
Sw = Vehicle_Data(23);
t = Vehicle_Data(24);


end