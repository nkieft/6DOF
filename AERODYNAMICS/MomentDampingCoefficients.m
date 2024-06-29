function [C_l_p, C_m_q, C_n_r] = MomentDampingCoefficients(C_N_V_D, C_N_W_D, beta, alpha, L_N, L_F, L_T, CG, R_B, Vel, P, Q, R, dmdt, A_ref, q)
%MOMENTDAMPINGCOEFFICIENTS calculates moments produced by localized angles
%of attack on components caused by vehicle movement, which tend to oppose
%sed motion. DOES NORMAL FORCE ONLY ATM, NOT ROLL DAMPING

%Convert AoA to radians
beta = beta*pi/180;
alpha = alpha*pi/180;

%Convert CG to m
CG = CG/39.37;

%Calculates CNA from CN and alpha
C_N_alpha_V_D = abs(C_N_V_D/beta);
C_N_alpha_W_D = abs(C_N_W_D/alpha);

%Replaces NaN (from stationary vehicle)
C_N_alpha_V_D(isnan(C_N_alpha_V_D)) = 0;
C_N_alpha_W_D(isnan(C_N_alpha_W_D)) = 0;

%Lever arms from fins and nosecone about CG
r_F = abs(CG-(L_T-L_F/2));
r_N = abs(CG-(2/3*L_N));

%Damping derivatives: from apogee and aspire space articles:
%https://apogeerockets.com/education/downloads/Newsletter195.pdf
%http://www.aspirespace.org.uk/downloads/A%20dynamic%20stability%20rocket%20simulator.pdf
C_m_q_a = -Q/Vel*(C_N_alpha_W_D(1)*r_N^2+C_N_alpha_W_D(2)*r_F^2)/(R_B*2);
C_m_r_a = -R/Vel*(C_N_alpha_V_D(1)*r_N^2+C_N_alpha_V_D(2)*r_F^2)/(R_B*2);

%SET TO ZERO FOR NOW ONG
%Jet Damping: Apogee source above
C_l_p_j = -P*abs(dmdt)*(R_B)^2/A_ref/R_B/2/q;
C_m_q_j = -R*abs(dmdt)*(L_T-CG)^2/A_ref/R_B/2/q;
C_n_r_j = -Q*abs(dmdt)*(L_T-CG)^2/A_ref/R_B/2/q;

%Net damping coefficients (sum of aerodynamic and jet)
C_l_p = C_l_p_j;
C_m_q = C_m_q_a+C_m_q_j;
C_n_r = C_m_r_a+C_n_r_j;

end