function ECEF = ECIToECEF(t, ECI, lon)
%ECITOECEF transforms vectors from the inertial to rotating earth reference
%frame. It does not accound for rotation, only transforms the vector.

%Earth rotation rate
omega = 7.2921159*10^-5;

%Angle rotated since launch
theta = t*omega-sign(t)*deg2rad(lon);
theta_deg = rad2deg(theta);
%Transformation matrix [the lords only use for APPM 3310]
T = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
T = T';

%Output vector
ECEF = T*ECI;

end