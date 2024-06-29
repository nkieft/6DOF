function [Body] = ECIToBody(ECI, q0, q1, q2, q3)
%Transforms a vector from the inertial frame to the body frame. Does not
%account for motion of the frames relative to one another. 

%Transformation matrix (from quats)
T = [2*(q0^2+q1^2)-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2); ...
    2*(q1*q2+q0*q3), 2*(q0^2+q2^2)-1, 2*(q2*q3-q0*q1); ...
    2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(q0^2+q3^2)-1];
T = transpose(T);

%Output vector
Body = T*ECI;

end