function init = FillInit(V_dat, L_dat)
%FILLINIT fills the initial state vector for the ascent function. See
%ASCENT for more details on state vector allocations.

%Body velocities
U0 = 0;
V0 = 0;
W0 = 0;

%Body rotations
P0 = 0;
Q0 = 0;
R0 = 0;

%Initial attitude
theta_0 = -deg2rad(L_dat(1));
psi_0 = deg2rad(L_dat(4));
phi_0 = pi/2-deg2rad(L_dat(3));
[q0_0, q1_0, q2_0, q3_0] = EulerToQuat(theta_0, phi_0, psi_0);

%Initial ECEF position
lat = L_dat(1);
lon = L_dat(2);
h = L_dat(8);
[Xe_0, Ye_0, Ze_0] = GeoToECEF(lat, lon, h);

%Initial ECEF velocity
Vxe_0 = 0;
Vye_0 = 0;
Vze_0 = 0;

%Initial ECI position
[Xi_0, Yi_0, Zi_0] = GeoToECEF(lat, 0, h);
%Initial ECI velocity
Vxi_0 = 0;
Vyi_0 = Xi_0*7.292115* 10^-5;
Vzi_0 = 0;

%Initial mass
m_0 = V_dat(4);

%Initial dist traveled
Uint_0 = 0;

%Put in a nice vector!
init = [U0; V0; W0; P0; Q0; R0; q0_0; q1_0; q2_0; q3_0; Xi_0; Yi_0; Zi_0; Vxi_0; Vyi_0; Vzi_0; m_0; Xe_0; Ye_0; Ze_0; Vxe_0; Vye_0; Vze_0; Uint_0];


end