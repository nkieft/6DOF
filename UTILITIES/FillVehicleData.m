function const = FillVehicleData(R_v, L_v, CG_0, CG_F, m_0, m_f, I_tot, C_r, L_n, ymac, cant, h_m, am, cdm, ad, cdd, I_i, I_f, T_a, T_o, C_t, Sw, t_F, b, X_F)
%FILLVEHICLEDATA consolidates useful vehicle parameters into a single
%vector! Also performs unit conversions.Structures were considered, but not 
%used due to how sloooooow they are.

    const(1) = (R_v/39.37)^2*pi; %Aerodynamic reference area [m^2]
    const(2) = CG_0; %Initial CG position [in]
    const(3) = CG_F; %Burnout CG position [in]
    const(4) = m_0/2.207; %Initial mass [kg]
    const(5) = m_f/2.207; %Final mass [kg]
    const(6) = I_tot; %Total motor impulse [N-s]
    const(7) = L_v/39.37; %Total vehicle length [m]
    const(8) = R_v/39.37; %Vehicle raidus [m]
    const(9) = C_r/39.37; %Fin length / root chord [m]
    const(10) = L_n/39.37; %Nosecone length
    const(11) = ymac/39.37; %Mean aerodynamic chord along lateral axis [m]
    const(12) = deg2rad(cant); %Fin cant angle [rad]
    const(13) = h_m; %Main deploy height [m]
    const(14) = am; %Main chute area [m^2]
    const(15) = cdm; %Main chute Cd [ul]
    const(16) = ad; %Drogue chute area [m^2]
    const(17) = cdd; %Drogue chute Cd [ul]
    const(18) = I_i; %Initial lateral inertia [kg-m^2]
    const(19) = I_f; %Final lateral inertia [kg-m^2]
    const(20) = T_a;
    const(21) = T_o;
    const(22) = C_t / 39.37;
    const(23) = Sw / 39.37;
    const(24) = t_F / 39.37;
    const(25) = b / 39.37;
    const(26) = X_F/39.37;

end