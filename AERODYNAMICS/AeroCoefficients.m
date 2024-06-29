function [C_Y, C_N, C_A, alpha_tot, beta, alpha_p] = AeroCoefficients(alpha_tot, alpha_p, beta, C)
%AEROCOEFFICIENTS uses Rasaero lookup tables to find aerodynamic
%coefficients along all body, and also calculates angle of attack (both net
%and about body y and z axes [deg]

%Side Force Coefficient Lookup
beta_abs = abs(beta); %For lookup
C_Y = AeroData(C(:, :, 1), beta_abs);

%Normal force coefficient lookup about body z axis
alpha_p_abs = abs(alpha_p); %For lookup
C_N= AeroData(C(:, :, 1), alpha_p_abs);

%Axial force coefficient (force along body x axis)
C_A= AeroData(C(:, :, 2), alpha_tot);

end