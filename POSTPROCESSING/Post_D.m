function output = Post_D(t, s, Vehicle_Data, Launchsite_Data, Physical_Const, W_h)
%POST_D is a modified descent file to get useful outputs

dSdt = zeros(6, 1); %Pre-allocate

X_i = s(1); %Inertial X position [m]
Y_i = s(2); %Inertial Y position [m]
Z_i = s(3); %Inertial Z position [m]
Vx_i = s(4); %Inertial X velocity [m/s]
Vy_i = s(5); %Inertial Y velocity [m/s]
Vz_i = s(6); %Inertial Z velocity [m/s]
X_e = s(7); %ECEF X position [m]
Y_e = s(8); %ECEF Y position [m]
Z_e = s(9); %ECEF Z position [m]
Vx_e = s(10); %ECEF X velocity [m/s]
Vy_e = s(11); %ECEF Y velocity [m/s]
Vz_e = s(12); %ECEF Z velocity [m/s]


omega = [0; 0; Physical_Const(1)]; %Earth rotation rate [rad/s]
H_m = Vehicle_Data(13)+Launchsite_Data(8); %Main deploy height [m]
Am = Vehicle_Data(14); %Main chute area [m^2]
Cdm = Vehicle_Data(15); %Main CD [ul]
Ad = Vehicle_Data(16); %Drogue chute area [m^2]
Cdd = Vehicle_Data(17); %Drogue CG [ul]
W_s = Launchsite_Data(6); %Wind speed [m/s]
W_az = Launchsite_Data(7); %Wind azimuth [deg]
m = Vehicle_Data(5); %Vehicle mass [kg]
L_long = Launchsite_Data(2);

%Calculate vehicle elevation
elevation = norm([X_i, Y_i, Z_i])-6378100;

%Choose which parachute to use based on elevation
[A, Cd] = SelectParachute(elevation, H_m, Ad, Cdd, Am, Cdm);

%Basic atmospheric/freestream properties
[T_a, a, P_a, rho_a] = atmosisa(elevation);

%Wind velocity vector in ECI frame
[X_w, Y_w, Z_w] = CalculateWind(W_s, W_az, X_i, Y_i, Z_i, X_e, Y_e, Z_e, elevation, omega, t, true, Launchsite_Data(1), Launchsite_Data(2), W_h);

%Difference between vehicle and atmosphere velocity in ECI frame
[Vx_eff, Vy_eff, Vz_eff] = EffectiveVelocityD(Vx_i, Vy_i, Vz_i, X_w, Y_w, Z_w);

%Speed of vehicle
speed = norm([Vx_eff, Vy_eff, Vz_eff]);

%Calculate dynamic pressure
q = DynamicPressure(rho_a, speed);

%Put velocity in vector
Vel = [Vx_eff; Vy_eff; Vz_eff];

%Calculate gravity force
F_G = -[X_i; Y_i; Z_i]/norm([X_i; Y_i; Z_i])*m*9.81;

%Calculate drag force
F_D = -A*Cd*q*Vel/speed;

%Sum forces, determine acceleration in ECI frame
F = F_G+F_D;
A_i = F/m;

%Correct for rotation of ECEF frame to find acceleration in it
A_e_i = ECIToECEF(t, A_i, L_long);
rotation_corr = - cross(2*omega, [Vx_e; Vy_e; Vz_e])-cross(omega, cross(omega, [X_e; Y_e; Z_e]));
A_e = A_e_i+rotation_corr;

[lon, lat, h] = ECEFtoGeo(X_e, Y_e, Z_e);

%Assign Output
output = zeros(1, 12);
output(1) = t;
output(2) = 0;
output(3) = 0;
output(4) = 0;
output(5) = 0;
output(6) = 0;
output(7) = norm(F_D);
output(8) = 0;
output(9) = speed/a;
output(10) = lat;
output(11) = lon;
output(12) = m;
output(13) = h*3.281;
output(14) = 1;
output(15) = 0;
output(16) = 0;
output(17) = 0;
output(18) = speed;
output(19) = norm(A_i);
output(20) = q;
output(21) = 0;
output(22) = 0;
output(23) = 0;
output(24) = 0;
output(25) = 0;
output(26) = 0;
output(27) = 0;
output(28) = 0;
output(29) = 0;
end