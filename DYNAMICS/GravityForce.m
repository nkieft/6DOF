function F_G = GravityForce(X_i, Y_i, Z_i, q0, q1, q2, q3, m, Uint, elevation)
%GRAVITYFORCE calulates a gravitational force on the vehicle in the direction
%of the center of the Earth. 

%Position unit vector vector in inertial frame
r = [X_i; Y_i; Z_i]/norm([X_i, Y_i, Z_i]);

%Position unit vector in body frame (direction of gravity)
r_b = ECIToBody(r, q0, q1, q2, q3);

%Variable gravity formula for weight force
W = 6.6743*10^-11*m*5.97219E24/(6378100+elevation)^2;

%Only apply gravity when rocket has started moving (to prevent falling thru
%ground)
if(Uint > 0.001)

    %Weight force towards earth CG
    F_G = -1*[W*r_b(1); W*r_b(2); W*r_b(3)];

else
    %Gravity reacted by ground/rail prior to launch
    F_G = [0;0;0];
end
end