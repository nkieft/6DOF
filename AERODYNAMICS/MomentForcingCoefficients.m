function [C_l, C_m, C_n] = MomentForcingCoefficients(C_Y_D, C_N_D, beta, alpha, F_cant, R_B, MAC_y, L_N, L_F, L_T, CG, P, U)
%MOMENTFORCINGCOEFFICIENTS calculates the moments that arise as a result of
%the vehicle being at an angle of attack w/ respect to free-stream. These
%moment coefficients are given in terms of a reference length equal to the
%diameter of the vehicle. 

%CG conversion to m
CG = CG/39.37;

%Angle of attacks around body y and z axes
beta = beta*pi/180;
alpha = alpha*pi/180;

%Calculates C_N_alpha for nosecone from distribution and AoA
C_Y_beta_N = abs(C_Y_D(1)/beta);
C_N_alpha_N = abs(C_N_D(1)/alpha);

%Calculates C_N_alpha for fins from distribution and AoA
C_Y_beta_F = abs(C_Y_D(2)/beta);
C_N_alpha_F = abs(C_N_D(2)/alpha);

%Accounts for NaNs (mostly from stationary vehicle)
if(isnan(C_Y_beta_F)), C_Y_beta_F = 0; end
if(isnan(C_N_alpha_F)), C_N_alpha_F = 0; end
if(isnan(C_Y_beta_N)), C_Y_beta_N = 0; end
if(isnan(C_N_alpha_N)), C_N_alpha_N = 0; end

%Mean CNA used for rolling moment calculations
C_N_alpha_F = (C_N_alpha_F+C_Y_beta_F)/4;

%Calculates effective fin cant (modified by roll rate)
F_cant_eff = F_cant-atan(P*(R_B+MAC_y)/U); %Effective cant
if(isnan(F_cant_eff)), F_cant_eff = F_cant; end %Remove nan

%Calculate roll moment coefficient from C_N_alpha
C_l = 4*C_N_alpha_F*F_cant_eff*(R_B+MAC_y)/R_B/2; 

%Lever arms of nosecone and fins about CG
LN = 2/3*L_N-CG;
LF = L_T-L_F/2-CG;

%Moment coefficient contributions from nosecone and fins, and sum about
%vehicle y axis
C_n_N = LN*C_Y_D(1)*sign(beta);
C_n_F = LF*C_Y_D(2)*sign(beta);
C_n = -(C_n_N+C_n_F)/R_B/2;

%Moment coefficient contributions from nosecone and fins, and sum about
%vehicle z axis
C_m_N = LN*C_N_D(1)*sign(alpha);
C_m_F = LF*C_N_D(2)*sign(alpha);
C_m = (C_m_N+C_m_F)/R_B/2;

end