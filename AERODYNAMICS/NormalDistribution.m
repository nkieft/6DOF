function [C_Y_D, C_N_D] = NormalDistribution(C_Y, C_N, CP, L_T, L_N, L_F, X_F)
%NORMALDISTRIBUTION estimates the distribution of the normal force
%coefficient between the fins and nosecone about both vehicle y and z axes.
%The first element of the output vectors is the nosecone C_N, and the
%second is the fins C_N

%Convert CP location to m
CP = CP/39.37;

%Matrix to solve: Solves system of equations that dictate that sum of
%normal force coefficients equals net (first row) and that the two forces
%will produce the desired center of pressure (second row)
A = [1, 1; 2/3*L_N-CP, L_T-(L_F-X_F)-CP];

%Right-hand-side matrices. Contain conditions established above
b_V = [C_Y; 0];
b_W = [C_N; 0];

%Solves equations for y and z axes to get force distrobutions.
C_Y_D = A\b_V;
C_N_D = A\b_W;

    

end