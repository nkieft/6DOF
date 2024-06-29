function [X_w, Y_w, Z_w] = CalculateWind(W_s, W_az, X_i, Y_i, Z_i, X_e, Y_e, Z_e, elevation, omega, t, inertial, L_lat, L_lon, W_h)
%CALCULATEWIND determines the velocity of the wind in terms of ECEF axes.
% See ASCENT for input and output variable descriptions

%Edge cases for extreme elevation sometimes encountered by runge kutta
if(elevation < 0)
    elevation = 0;
end
if(elevation > 500000)
    elevation = 500000;
end

h_res = W_h(2, 1)-1;

%Uses wind model (if selected by user) to determine wind in ENU frame
W_E = FastInterp(W_h(:, 1), W_h(:, 3), elevation, h_res);
W_N = FastInterp(W_h(:, 1), W_h(:, 2), elevation, h_res);
W_U = 0;

%Extract earth rotation rate
omega1 = omega(3);

%Wind vector in ENU coordinates
W_ENU = [W_E; W_N; W_U];

%Transforms ENU wind to ECEF frame
W_o = ENUToECEF(X_e, Y_e, Z_e, W_ENU);

%Inertial Conversion
if(inertial)
    W_o = ECIToECEF(-t, W_o, L_lon);
    W_o = W_o+cross(omega, [X_i; Y_i; Z_i;]);
end

X_w = W_o(1);
Y_w = W_o(2);
Z_w = W_o(3);

end