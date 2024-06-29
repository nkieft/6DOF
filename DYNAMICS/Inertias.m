function [I_A, I_N] = Inertias(m, r_v, I_0, I_f, m_0, m_f)
%INERTIAS calculates the lateral and axial MMOIs for the vehicle, assuming
%that it is a solid cylinder longitudinally. A correction factor is added. 

%%OLD
%Lateral MMOI, correction 0.6 added emperically
%I_N = (m/12)*(3*r_v^2+L_v^2)*0.6;

%%NEW, linear interpolation based on mass
I_N = I_f+(I_0-I_f)/(m_0-m_f)*(m-m_f);

%Axial MMOI, cylinder
I_A = m*r_v^2/2*1.5;

end