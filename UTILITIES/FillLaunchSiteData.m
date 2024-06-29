function const = FillLaunchSiteData(L_lat, L_long, L_Az, L_ang, R_L, L_e, RFC, MU)
%FILLLAUNCHSITEDATA consolidates data associated with launch conditions
%into a vector, and performs some unit conversions. 

const(1) = L_lat; %Launchsite latitude [deg]
const(2) = L_long; %Launchsite longitude [deg]
const(3) = L_Az; %Launcher azimuth [deg]
const(4) = L_ang; %Launcher angle [deg]
const(5) = R_L/3.281; %Rail length [m]
const(8) = L_e; %Launcher elevation [m]
const(10) = RFC*4.448; %Rail contact force [N]
const(11) = MU; %Rail friction coeff

end