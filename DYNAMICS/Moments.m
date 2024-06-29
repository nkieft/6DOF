function M = Moments(C_l, C_m, C_n, C_l_p, C_m_q, C_n_r, q, A_ref, L_ref, onrail, F_T, L_v, CG)
%MOMENTS combines damping and forcing coefficients to determine net moments
%on the vehicle

%Off-Axis Thrust Moment
r = L_v-CG/39.37;
C_m_T = F_T(2)*r/L_ref/A_ref/q;
C_n_T = F_T(3)*r/L_ref/A_ref/q;

%Combine damping and forcing moment coefficients, along with dynamic
%pressure, reference area, reference length to get moments about all body
%axes
M(1) = (C_l+C_l_p)*q*A_ref*L_ref;
M(2) = -(C_m-C_m_q+C_m_T)*q*A_ref*L_ref;
M(3) = -(C_n-C_n_r+C_n_T)*q*A_ref*L_ref;


end