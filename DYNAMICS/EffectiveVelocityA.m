function [Ue, Ve, We] = EffectiveVelocityA(Vx_e, Vy_e, Vz_e, X_w, Y_w, Z_w, t, q0, q1, q2, q3, lon)
%EFFECTIVEVELOCITYA Calculates the vehicle's effective velocity within the
%atmosphere along each of the body axes

%Calculates difference between current velocity and atmospheric velocity in
%ECEF frame
X_eff = Vx_e-X_w;
Y_eff = Vy_e-Y_w;
Z_eff = Vz_e-Z_w;

%Put in vector
vec = [X_eff; Y_eff; Z_eff];

%Transform from ECEF to ECI frames (using negative time :D)
ECI = ECIToECEF(-t, vec, lon);

%Transform from ECI frame to Body frame
UVW = ECIToBody(ECI, q0, q1, q2, q3);

%Assign individual components in X, Y, Z body directions
Ue = UVW(1);
Ve = UVW(2);
We = UVW(3);

end