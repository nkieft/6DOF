function [q0, q1, q2, q3] = EulerToQuat(phi, theta, psi)
%EULER2QUAT calculates the quaternions associated with a given series of
%euler angle rotations. Used for initial orientation of the vehicle only

%Idk why i kept this here but it'll stay i guess
%q1 = cos(psi/2)*cos(theta/2)*cos(phi/2)+sin(psi/2)*sin(theta/2)*sin(phi/2);
%q2 = sin(psi/2)*cos(theta/2)*cos(phi/2)-cos(psi/2)*sin(theta/2)*sin(phi/2);
%q3 = cos(psi/2)*sin(theta/2)*cos(phi/2)+sin(psi/2)*cos(theta/2)*sin(phi/2);
%q4 = cos(psi/2)*cos(theta/2)*sin(phi/2)-sin(psi/2)*sin(theta/2)*cos(phi/2);

%Quaterions (find source for math (i prolly read in a rando paper))
q0 = (sin(phi/2)*sin(theta/2)*sin(psi/2)+cos(phi/2)*cos(theta/2)*cos(psi/2));
q1 = (sin(phi/2)*cos(theta/2)*sin(psi/2)+cos(phi/2)*sin(theta/2)*cos(psi/2));
q2 = (sin(phi/2)*cos(theta/2)*cos(psi/2)-cos(phi/2)*sin(theta/2)*sin(psi/2));
q3 = (-sin(phi/2)*sin(theta/2)*cos(psi/2)+cos(phi/2)*cos(theta/2)*sin(psi/2));

end