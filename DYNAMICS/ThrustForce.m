function F_T = ThrustForce(t, T, F_T_a, F_T_o)
%THRUSTFORCE calculates the thrust produced by the motor [N] (from the
%supplied file) via linear interpolation with the current simulation time
%step

%If vehicle is in thrust phase
if(t < T(end, 1) && t > 0)

    F_T = FastInterp(T(:, 1), T(:, 2), t, 0.03);
    
else %Thrust 0 when burnout has happened
    F_T = 0;
end

if(isnan(F_T)) %usually happens after burnout
        F_T = 0;
end

F_T_x = F_T*cosd(F_T_a);
F_T_y = F_T*sind(F_T_a)*cosd(F_T_o);
F_T_z = F_T*sind(F_T_a)*sind(F_T_o);

%Put it in a vector
F_T = [F_T_x; F_T_y; F_T_z];

end