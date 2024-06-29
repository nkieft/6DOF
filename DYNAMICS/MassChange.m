function dmdt = MassChange(F_T, I_tot, m_0, m_f)
%DMDT calculates the rate of change of the vehicle's mass as a function of
%it's current thrust, and initial and final mass

m_p = m_0-m_f; %Propellant mass

%Rate of change of vehicle mass is proportional to current thrust over
%total impulse (can be derived with chain rule)
dmdt = -F_T(1)/I_tot*m_p; 


end