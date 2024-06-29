function [ECI] = BodyToECI(Body, q0, q1, q2, q3)
%BODYTOINERTIAL transforms a vector from the body frame to the inertial
%frame. It does not accound for rotation of the frames. 

%Transformation matrix from quaternions 
T = [2*(q0^2+q1^2)-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2); ...
    2*(q1*q2+q0*q3), 2*(q0^2+q2^2)-1, 2*(q2*q3-q0*q1); ...
    2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(q0^2+q3^2)-1];

%Output vector
ECI = T*Body;

end