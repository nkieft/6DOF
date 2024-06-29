function CG = CenterOfGravity(m, CG_0, CG_f, m_0, m_f)
%CENTEROFGRAVITY calculates the vehicle's current CG distance from the tip
%[in] by linearly interpolating with it's current mass

%Linearly interpolate CG with mass (not perfect but sufficient)
CG = CG_f+(CG_0-CG_f)/(m_0-m_f)*(m-m_f);

end