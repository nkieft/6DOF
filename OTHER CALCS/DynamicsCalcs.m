function [C1, C2, omega_n, zeta, omega_n_c, zeta_c, AR_c_res, P, t] = DynamicsCalcs(hM, speed, V_D, I_N, I_A, C_N, C_N_N, C_N_F, alpha, x_cg, x_cp, t, P)

%% Derivations available in "Topics in Advanced Model Rocketry"

%Organize Geometric Definitions
x_cg = x_cg/39.37;
x_cp = x_cp/39.37;
L_T = V_D(7);
L_F = V_D(9);
L_N = V_D(10);
x_F = L_T-L_F/3;
x_N = 2/3*L_N;
A_ref = V_D(1);
alpha = deg2rad(alpha);

%Atmospheric Lookup
[~, ~, ~, rho] = atmosisa(hM, extended=true);

%Determining aero AoA derivatives
C_N_alpha = sign(alpha).*C_N./alpha;
C_N_N_alpha = sign(alpha).*C_N_N./alpha;
C_N_F_alpha = sign(alpha).*C_N_F./alpha;

%Forcing and damping coefficients
C1 = 1/2*rho.*speed.^2.*A_ref.*C_N_alpha.*(x_cp-x_cg);
C2 = 1/2*rho.*speed.*A_ref.*(C_N_N_alpha.*(x_N-x_cg).^2+C_N_F_alpha.*(x_F-x_cg).^2);

%Vehicle natural frequency, damping ratio, damped frequency
omega_n = sqrt(C1./I_N);
zeta = C2./(2*sqrt(C1.*I_N));
omega_d = omega_n.*sqrt(1-zeta.^2);

%Coupled ratios (No physical Meaning)
omega_n_c = sqrt(C1./(I_N+I_A));
zeta_c = C2./(2.*sqrt(C1.*(I_N+I_A)));
beta = 0:0.01:2;
AR_c_res = 1./(2.*zeta_c.*sqrt(1-zeta_c.^2));

end