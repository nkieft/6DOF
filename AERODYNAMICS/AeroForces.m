function F_A = AeroForce(A_ref, q, C_A, C_Y, C_N, F_T, onRail, dir, elev, elev_i)
%AEROFORCE calculates the aerodynamic forces on the vehicle from provided
%coefficients and freestream conditions

%Ensure vehicle is moving
if(q > 0 && elev>elev_i)

    %Calculate forces in each direction
    F_x = AeroForce(A_ref, q, C_A, -dir(1));
    F_y = -AeroForce(A_ref, q, C_Y, dir(2));
    F_z = -AeroForce(A_ref, q, C_N, dir(3));

%Forces zero when vehicle stationary
else
    F_x = 0;
    F_y = 0;
    F_z = 0;
end

%Rail reacts normal forces when rocket is on it
if(onRail)
    F_y = 0;
    F_z = 0;
else

end

%Output vector
F_A = [F_x; F_y; F_z];


end