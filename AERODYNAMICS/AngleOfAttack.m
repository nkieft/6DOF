function [alpha_tot, beta, alpha_p] = AngleOfAttack(U, V, W, onRail)

%Assume lateral velocity is zero when on the rail
if(onRail)
    V = 0;
    W = 0;
end

%Angle of attack about body yaw axis
beta = atand(V/U);
if(isnan(beta)), beta = 0; end %For stationary vehicle

%Angle of attack about body pitch axis
alpha_p = atand(W/U);
if(isnan(alpha_p)), alpha_p = 0; end %For stationary vehicle

%Net angle of attack of vehicle
alpha_tot = abs(acosd(dot([U, V, W], [U, 0, 0])/(norm([U, V, W])*U)));
if(isnan(alpha_tot)), alpha_tot = 0; end %For stationary vehicle


end