function F_A = AeroForce(A_ref, q, C_A, dir)
%DRAGFORCE applies the drag formula to find axial force on the vehicle
    
    %Calculate it yo
    F_A = dir*C_A*q*A_ref;
    
end