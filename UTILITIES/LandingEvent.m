function [value, isterminal, direction] = LandingEvent(t, s, L_data)
%LANDINGEVENT is the ode45 termination for the descent function. Triggers
%when vehicle goes beneath the ground

%Elevation of vehicle
[~, ~, elevation] = ECEFtoGeo(s(7), s(8), s(9));

%Did it land
value = elevation < L_data(8);

isterminal = 1;
direction = 0;

end