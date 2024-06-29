function F_R = RailForces(F_G, F_A, F_T, U, onRail, mu, F_c, m, l_a)

    if(onRail)
        F_R(2) = -F_G(2)-F_A(2)-F_T(2);
        F_R(3) = -F_G(3)-F_A(3)-F_T(3);
        F_c = F_c+m*sind(l_a);
        if(U > 0)
            F_R(1) = -mu*F_c;
        end
        F_R = F_R';

    else
        F_R = zeros(3,1);
    end
end