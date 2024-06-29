function [X_eff, Y_eff, Z_eff] = EffectiveVelocityD(Vx, Vy, Vz, X_w, Y_w, Z_w)
%EFFECTIVEVELOCITYA Calculates the vehicle's effective velocity within the
%atmosphere in ECEF frame

%Take differences
X_eff = Vx-X_w;
Y_eff = Vy-Y_w;
Z_eff = Vz-Z_w;

end