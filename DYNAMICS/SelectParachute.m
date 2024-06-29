function [A, CD] = SelectParachute(elevation, H_m, ad, cdd, am, cdm)
%SELECTPARACHUTE determines which chute area and cd to use, based on when
%the main should deploy and the current altitude of the vehicle.

%Use main below main deploy height
if(H_m>elevation)
    A = am;
    CD = cdm;

%Use drogue above main deploy height
else
    A = ad;
    CD = cdd;
end
end
