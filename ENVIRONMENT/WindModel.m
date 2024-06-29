function W = WindModel(L_lat, L_long, L_day, AV_Wind, AV_Wind_d, ManW, MaxA)
%WINDMODEL calculates wind data from sea level to 100 km in 10m increments
%
%Inputs: Launcher Latitude, Launcher Longitude, Launch Day
%
%Outputs: Wind speeds in north and east directions

%Heights from 1m to 100km in 50 m increments
res = 50;
h = 1:res:MaxA;

%Empty array for other inputs
temp = ones(length(h), 1);

%Wind data for all altitudes at launchsite
W = atmoshwm(temp*L_lat, temp*L_long, h, 'day', temp*L_day);
W(:, 3) = W(:, 2);
W(:, 2) = W(:, 1);
W(:, 1) = h;

if(ManW)
alt_AV = [0 3 6 9 12 15 18 21 24 27 30 36 42 48]/3.281*1000;

alt_AV_q = 0:res:alt_AV(end);

AV_w_speed = interp1(alt_AV, AV_Wind, alt_AV_q)/1.944;
AV_w_dir = interp1(alt_AV, AV_Wind_d, alt_AV_q);

AV_N = -AV_w_speed.*cosd(AV_w_dir);
AV_E = -AV_w_speed.*sind(AV_w_dir);

W(1:length(AV_N), 2) = AV_N;
W(1:length(AV_E), 3) = AV_E;
end

end