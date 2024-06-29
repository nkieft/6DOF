function CP = CenterOfPressure(alpha_tot, Cp_data, q)
%Calculates the vehicle's center of pressure [in] distance from the tip
%based on vehicle AoA and mach number. (Data from RASAero)

%Calculates angle between velocity and orientation vectors
AoA = alpha_tot;

%Mostly accounts for when vehicle is still stationary
if(isnan(AoA)), AoA = 0; end

%Only vary CP when the vehicle is going sorta fast. Otherwise, use initial
%value
%This is pretty much only for when it's on the rail
if(q > 1000)
    CP = AeroData(Cp_data, AoA);
else
    CP = Cp_data(1, 1);
end

end