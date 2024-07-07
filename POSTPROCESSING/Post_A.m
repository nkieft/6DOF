function output = Post_A(t, s, V_data, L_data, Const, F_T_data, C_b_data, C_c_data, Cp_data, W_h, atm_data)
%POST_A is a modified ascent file to get useful outputs

%% Extract State Variables
%----------------------------------------------------------------------
U = s(1); %Velocity along body x axis [m/s]
V = s(2); %Velocity along body y axis [m/s]
W = s(3); %Velocity along body z axis [m/s]
P = s(4); %Angular velocity about body x axis [rad/s]
Q = s(5); %Angular velocity about body y axis [rad/s]
R = s(6); %Angular velocity about body z axis [rad/s]
q0 = s(7); %Quaterion 0 
q1 = s(8); %Quaterion 1
q2 = s(9); %Quaterion 2
q3 = s(10); %Quaterion 3
X_i = s(11); %Inertial X position [m]
Y_i = s(12); %Inertial Y position [m]
Z_i = s(13); %Inertial Z position [m]
Vx_i = s(14); %Inertial X velocity [m/s]
Vy_i = s(15); %Inertial Y velocity [m/s]
Vz_i = s(16); %Inertial Z velocity [m/s]
m = s(17); %Vehicle mass [kg]
X_e = s(18); %ECEF X position [m]
Y_e = s(19); %ECEF Y position [m]
Z_e = s(20); %ECEF Z position [m]
Vx_e = s(21); %ECEF X velocity [m/s]
Vy_e = s(22); %ECEF Y velocity [m/s]
Vz_e = s(23); %ECEF Z velocity [m/s]
Uint = s(24); %Distance travelled along body x axis [m]

%----------------------------------------------------------------------

%% Extract Constants
%----------------------------------------------------------------------
L_r = L_data(5); %Rail length [m]
A_ref = V_data(1); %Reference area [m^2]
CG_0 = V_data(2); %CG distance from tip initially [in]
CG_f = V_data(3); %CG distance from tip at burnout [in]
m_0 = V_data(4); %Initial vehicle mass [kg]
m_f = V_data(5); %Burnout vehicle mass [kg]
I_tot = V_data(6); %Total motor impulse [N-s]
r_v = V_data(8); %Vehicle radius [m]
L_v = V_data(7); %Vehicle length [m]
omega = [0; 0; Const(1)]; %Earth rotation rate vector [rad/s]
W_s = L_data(6); %Wind speed [m/s]
W_az = L_data(7); %Wind azimuth [deg]
C_r = V_data(9); %Fin root chord [m]
L_N = V_data(10); %Nosecone length [m]
ymac = V_data(11); %Fin aerodynamic chord [m]
F_cant = V_data(12); %Fin cant angle [rad]
L_lat = L_data(1); %Launcher latitude [deg]
L_long = L_data(2); %Launcher longitude [deg]
I_0 = V_data(18); %Initial lateral inertia [kg-m^2]
I_f = V_data(19); %Final lateral inertia [kg-m^2]
h_0 = L_data(8); %Launcher elevation;
F_T_a = V_data(20);
F_T_o = V_data(21);
F_c = L_data(10);
mu = L_data(11);
alpha_L = L_data(4);
X_f = V_data(26);
%-------------------------------------------------------------------------

%% Basic State Info
%----------------------------------------------------------------------
%If vehicle has thrust plume, use motor burning coefficient, otherwise, use
%coasting coefficients.
if(t < F_T_data(end, 1)) C_data = C_b_data; else C_data=C_c_data; end

%Determine if vehicle is on launch rail (which constrains it's motion)
onRail = Uint < L_r;

%Vehicle elevation
[~, ~, elevation] = ECEFtoGeo(X_e, Y_e, Z_e);
if(elevation < 0)
    pause(0.01);
    elevation = 0;
end
%----------------------------------------------------------------------

%% Freestream Properties
%----------------------------------------------------------------------
%Basic atmospheric properties [K, m/s, Pa, kg/m^3]
a = FastInterp(atm_data(:, 1), atm_data(:, 2), elevation, atm_data(2, 1)-1);
rho_a = FastInterp(atm_data(:, 1), atm_data(:, 3), elevation, atm_data(2, 1)-1);

%Calculate wind velocity in body frame [m/s]
[W_x, W_y, W_z] = CalculateWind(W_s, W_az, X_i, Y_i, Z_i, X_e, Y_e, Z_e, elevation, omega, t, false, L_lat, L_long, W_h);

%Determine velocity in body frame relative to wind/atmosphere [m/s]
[U_eff, V_eff, W_eff] = EffectiveVelocityA(Vx_e, Vy_e, Vz_e, W_x, W_y, W_z, t, q0, q1, q2, q3, L_long);

%Magnitude of vehicle velocity relative to wind [m/s]
speed = norm([U_eff, V_eff, W_eff]);

%Vehicle direction relative to wind [ul]
dir = [sign(U_eff), sign(V_eff), sign(W_eff)];

%Dynamic pressure experienced by vehicle [Pa]
q = DynamicPressure(rho_a, speed);

%Mach number of vehicle [ul]
Mach = speed/a;
%----------------------------------------------------------------------

%% Aerodynamic Coefficients
%----------------------------------------------------------------------

%Determine Angle of Attack along each axis
[alpha_tot, beta, alpha_p] = AngleOfAttack(U_eff, V_eff, W_eff, onRail);

[C_index] = MachIndex(Mach, size(C_b_data, 1));

%Vehicle center of pressure (relative to tip) location [m]
CP = CenterOfPressure(alpha_tot, Cp_data(C_index, :), q);

%Aerodynamic coefficients along all body axes [ul], net angle of attack
%[deg], and angles of attack about body y and z axes [deg]
%[C_Y, C_N, C_A] = AeroCoefficients(alpha_tot, alpha_p, beta, C_data(C_index, :, :));

%Yaw normal coefficient
C_Y = AeroData(C_data(C_index, :, 1), abs(beta));

%Pitch normal coefficient
C_N = AeroData(C_data(C_index, :, 1), abs(alpha_p));

%Axial force coefficient
C_A = AeroData(C_data(C_index, :, 2), alpha_tot);

%Normal force coefficients of vehicle along x and y axes, distributed among
%nosecone and fin sections [ul]
[C_Y_D, C_N_D] = NormalDistribution(C_Y, C_N, CP, L_v, L_N, C_r, X_f);
%----------------------------------------------------------------------

%% Calculate Forces
%----------------------------------------------------------------------
%Thrust force [N]
F_T = ThrustForce(t, F_T_data, F_T_a, F_T_o);

%Aerodynamic forces [N]
F_A = AeroForces(A_ref, q, C_A, C_Y, C_N, F_T, onRail, dir, elevation, h_0);

%Gravitational force [N]
F_G = GravityForce(X_i, Y_i, Z_i, q0, q1, q2, q3, m, Uint, elevation);

%Rail Reaction Forces
F_R = RailForces(F_G, F_A, F_T, U, onRail, mu, F_c, m, alpha_L);

%Net force is sum of above [N]
F = F_T+F_A+F_G+F_R;
%----------------------------------------------------------------------

%% Calculate Mass Rate of Change
%----------------------------------------------------------------------
%Mass rate of change calculated from thrust [kg/s]
dmdt = MassChange(F_T, I_tot, m_0, m_f);
%----------------------------------------------------------------------

%% Calculate Moments
%----------------------------------------------------------------------
%Current vehicle CG (relative to tip) [m]
CG = CenterOfGravity(m, CG_0, CG_f, m_0, m_f);

%Aero Restoring / Forcing moment coefficients about all body axes [ul]
[C_l, C_n, C_m] = MomentForcingCoefficients(C_Y_D, C_N_D, beta, alpha_p, F_cant, r_v, ymac, L_N, C_r, L_v, CG, P, U);

%Damping moment coefficients about body y and z axes [ul]
[C_l_p, C_m_q, C_n_r] = MomentDampingCoefficients(C_Y_D, C_N_D, beta, alpha_p, L_N, C_r, L_v, CG, r_v, speed, P, Q, R, dmdt, A_ref, q);

%Net moments about all vehicle axes [N-m]
M_F = Moments(C_l, C_n, C_m, C_l_p, C_m_q, C_n_r, q, A_ref, r_v*2, onRail, F_T, L_v, CG);

%Rail Reaction Moments
M_R = RailMoments(onRail, M_F);

M = M_F+M_R;
%----------------------------------------------------------------------

%% Calculate Intertias
%----------------------------------------------------------------------
%Moment of inertia about longitudinal and lateral axes [kg-m^2]
[I_A, I_N] = Inertias(m, r_v, I_0, I_f, m_0, m_f);
%----------------------------------------------------------------------

%% Calculate Derivatives
%----------------------------------------------------------------------
%Rotational accelerations about body x, y, z axes [rad/s/s]
P_dot = (M(1)-(I_N-I_N)*Q*R)/I_A;
Q_dot = (M(2)-(I_A-I_N)*R*P)/I_N;
R_dot = (M(3)-(I_N-I_A)*P*Q)/I_N;

%Translational accelerations in body x, y, z, axes [m/s/s]
U_dot = F(1)/m+R*V-Q*W;
V_dot = F(2)/m+P*W-R*U;
W_dot = F(3)/m+Q*U-P*V;

%Attitude rates of change in quaternions [1/s]
lambda = 1*(1-(q0^2+q1^2+q2^2+q3^2));
q0_dot = (q1*P+q2*Q+q3*R)/-2 + lambda*q0;
q1_dot = (q0*P+q2*R-q3*Q)/2 + lambda*q1;
q2_dot = (q0*Q+q3*P-q1*R)/2 + lambda*q2;
q3_dot = (q0*R+q1*Q-q2*P)/2 + lambda*q3;

%Acceleration vector in inertial frame [m/s/s]
A_i = BodyToECI(F/m, q0, q1, q2, q3);

%Transform inertial acceleration vector to ECEF 
A_e_i = ECIToECEF(t, A_i, L_long);

%Coerrection factor for rotation of earch
rotation_corr = - cross(2*omega, [Vx_e; Vy_e; Vz_e])-cross(omega, cross(omega, [X_e; Y_e; Z_e]));

%Acceleration vector in ECEF frame [m/s/s]
A_e = A_e_i+rotation_corr;

%----------------------------------------------------------------------

%Removes weird high AoA initially
if(alpha_tot > 90)
    alpha_tot = 180-alpha_tot;
end

%% Assign Output
%----------------------------------------------------------------------

[lon, lat, h] = ECEFtoGeo(X_e, Y_e, Z_e);

output = zeros(1, 12);
output(1) = t;
output(2) = alpha_tot;
output(3) = CG;
output(4) = CP;
output(5) = (CP-CG)/39.37/(2*r_v);
output(6) = sqrt(F(2)^2+F(3)^2);
output(7) = F_A(1);
output(8) = F_T(1);
output(9) = speed/a;
output(10) = lat;
output(11) = lon;
output(12) = m;
output(13) = h*3.281;
output(14) = q0;
output(15) = q1;
output(16) = q2;
output(17) = q3;
output(18) = speed;
output(19) = norm(A_i);
output(20) = q;
output(21) = P;
output(22) = Q;
output(23) = R;
output(24) = I_N;
output(25) = I_A;
output(26) = C_N;
output(27) = C_N_D(1);
output(28) = C_N_D(2);
output(29) = alpha_p;


%----------------------------------------------------------------------



end