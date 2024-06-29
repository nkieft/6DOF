function [x_cg, I, m_tot, dmdx] = MassDistribution(s, dx, res, x, MO, IO, CGO)
%MASSDISTRIBUTION Calculates the total mass, center of gravity, and
%Inertia of a vehicle configuration. 
%Author: Nick Kieft
%Contributor: Alexis Villarose

%Mass Distribution Array
dmdx = zeros(res, 1);

%% Mass Distribution
for i = 1:height(s)

%Find index of relevant elements
L = s.XF(i)-s.X0(i);
index_start = round(s.X0(i)/dx)+1;
index_end = round(s.XF(i)/dx);
elements = index_end-index_start+1;

%Find radii of relevant elements
R0 = s.D0(i)/2;
RF = s.DF(i)/2;

%Area of trapezoidal dist. and rectangular contribution
a_tp = (R0+RF)/2*L;
a_r = R0*L;

%Find component contribution to mass distribution
m = s.M(i);
m_r = a_r/a_tp*m;
h_r = m_r/L;
m_t = m-m_r;
h_t = 2*m_t/L;
const = ((h_t+h_r)-h_r)/L;

%Create vector of mass distribution of component
x_temp = linspace(0, L, elements)';
dmdx_s = (h_r+const.*x_temp);

%Add mass distribution to vector
dmdx(index_start:index_end) = dmdx(index_start:index_end)+dmdx_s;

end

%% CG Integral
x_cg = trapz(x, dmdx.*x')./trapz(x, dmdx);

%% Inertia Integral
r = x-x_cg;
I = trapz(x, dmdx.*r'.^2);

%% Total Mass
m_tot = trapz(x, dmdx);

%Conversions
x_cg = x_cg*39.37;
m_tot = m_tot*2.207;

%% Overrides
if(CGO > 0)
    x_cg = CGO*39.37;
end

if(MO > 0)
    m_tot = MO*2.207;
end

if(IO > 0)
    I = IO;
end

end