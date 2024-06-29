function [V_D, V_M, QTYPins, W_BP, F_sep] = RecoveryCalcs(Vehicle_Data, Launchsite_Data, dmdx, x_c, Mach, q, l_c, m)

%Vehicle Data Extraction
L_T = Vehicle_Data(7);
L_F = Vehicle_Data(8);
L_N = Vehicle_Data(10);
L_B = L_T-L_F-L_N;
S_ref = Vehicle_Data(1);

x = linspace(0, L_T, 1000);
dx = x(2);

%Extraction of RASAero Data
data = readmatrix("CDDataFile0.txt", "Delimiter", ["   ", "  "]);
data(isnan(data(:, 1)), :) = [];
data_A = data(1:6:end, :);
C_Prot = data_A(1, 11);
M_index = round(Mach*100);
data_A = data_A(M_index, :);

%Index Calculation
i_c = round((x_c/39.377) / dx);
i_N = round(L_N / dx);
i_F = round(L_T / dx);
i_B = round((L_B+L_N) / dx);

% Axial Distribution
C_A_tot = data_A(3);
F_A_Tot = -C_A_tot*S_ref*q;
W_A = zeros(size(x));

%Nosecone Coeff
C_A_N = L_N/L_T*data_A(5)+data_A(6);
const = C_A_N*2/L_N^2;
W_A(1:i_N) = W_A(1:i_N)+const*x(1:i_N);

%Body Coeff
C_A_B = L_B/L_T*data_A(5)+data_A(12);
W_A(i_N+1:i_B) = W_A(i_N+1:i_B)+C_A_B/L_B;

%Fins
C_A_F = L_F/L_T*data_A(5)+data_A(8)+data_A(9)+data_A(10)+data_A(11);
W_A(i_B+1:i_F) = W_A(i_B+1:i_F)+C_A_F/L_F;

%Base Drag
C_Base = data_A(7);

m_u = trapz(x(1:i_c), dmdx(1:i_c));
m_l = trapz(x(1+i_c:end), dmdx(1+i_c:end));

F_u = trapz(x(1:i_c), W_A(1:i_c))*S_ref*q;
F_l = (trapz(x(1+i_c:end), W_A(1+i_c:end))+C_Base+C_Prot)*S_ref*q;

a_cg = abs(F_A_Tot/(m_u+m_l));
F_sep = F_l-(F_l+F_u)*m_l/(m_l+m_u);

QTYPins = ceil(F_sep/4.448/40*1.5);

V = l_c*S_ref;
F = QTYPins*40;
P = F/S_ref;
W_BP = 3*454*P*V/266/3307;

%Chute Data Extraction
h_m = Vehicle_Data(13);
A_m = Vehicle_Data(14);
Cd_m = Vehicle_Data(15);
h0 = Launchsite_Data(8);
A_d = Vehicle_Data(16);
Cd_d = Vehicle_Data(17);

[~, ~, ~, rho_D] = atmosisa(h_m+h0);
[~, ~, ~, rho_L] = atmosisa(h0);

V_D = sqrt(2*m*9.81/rho_D/A_d/Cd_d);
V_M = sqrt(2*m*9.81/rho_L/A_m/Cd_m);

end

