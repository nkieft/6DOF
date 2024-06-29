function [value, isterminal, direction] = ApogeeEvent(t, s, L_data, Ballistic)
%APOGEEEVENT is an ode45 event function. For ballistic recovery, the
%simulation terminates when the vehicle crosses the earth's radius. For
%chuted recovery, the simulation terminates when the vehicle has a significant
% velocity component in the direction towards the earth.

%Vehile elevation
[~, ~, elevation] = ECEFtoGeo(s(18), s(19), s(20));


if(Ballistic)
    %Vehicle back on earth
    value = elevation < L_data(8)-5;

else
    %Vehicle at apogee
    value = dot([s(18), s(19), s(20)], [s(21), s(22), s(23)]) < -5;
end

if(value == 1)
    pause(0.0001);
end

isterminal = 1;
direction = 0;

end