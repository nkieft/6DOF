function [W_A, W_N, F_A_Tot, F_N_Tot] = AeroDistribution(V_dat, Mach, q, F_T, filename)

A_ref = V_dat(1);
L_N = V_dat(10);
L_F = V_dat(9);
L_T = V_dat(7);
L_B = L_T-L_F-L_N;

x = linspace(0, L_T, 1000)';
dx = x(2);
res = length(x);

%% Data Formatting
data = readmatrix(filename, "Delimiter", ["   ", "  "] );
data(isnan(data(:, 1)), :) = [];

data_A = data(1:6:end, :);
data_N = data(4:6:end, :);

index = round(Mach*100);
data_A = data_A(index, :);
data_N = data_N(index, :);

%% Indexing
W_N = zeros(res, 1);
W_A = zeros(res, 1);

i_N = round(L_N/dx);
i_B = round((L_B+L_N)/dx);
i_F = round((L_N+L_B+L_F)/dx);

%% Axial Distribution
C_A_Tot = data_A(3);
F_A_Tot = C_A_Tot*A_ref*q;
F_A_Tot = -F_A_Tot+F_T;

%Nosecone
C_A_N = L_N/(L_T)*data_A(5)+data_A(6);
const = C_A_N*2/L_N^2;
W_A(1:i_N) = W_A(1:i_N)+const*x(1:i_N);

%Body
C_A_B = data_A(5)*(L_B/L_T)+data_A(12);
W_A(i_N+1:i_B) = W_A(i_N+1:i_B)+C_A_B/L_B;

%Fins
C_A_F = L_F/L_T*data_A(5)+data_A(8)+data_A(9)+data_A(10)+data_A(11);
W_A(i_B+1:i_F) = W_A(i_B+1:i_F)+C_A_F/L_F;

%% Normal Distribution
C_N_Tot = data_N(3);
F_N_Tot = C_N_Tot*A_ref*q;
CP = data_N(4)/39.37;

%Distribution
A = [1, 1; 2/3*L_N-CP, L_N+L_B+L_F/2-CP];
b = [C_N_Tot; 0];
c = A\b;


%Nosecone
const = 2*c(1)/L_N^2;
W_N(1:i_N) = W_N(1:i_N)+const*x(1:i_N);

%Fins
W_N(i_B:i_F)=W_N(i_B:i_F)+c(2)/L_F;

W_A = -W_A*A_ref*q;
W_N = W_N*A_ref*q;

end
