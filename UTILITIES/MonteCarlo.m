function [F_T, Mass_Initial, Mass_Final, Launcher_Azimuth, Launcher_Angle, C_b, C_c, F_Cant, F_T_a, F_T_o] = MonteCarlo(F_T, Mass_Initial, Mass_Final, Launcher_Azimuth, Launcher_Angle, C_b, C_c, F_Cant, F_T_a, F_T_o, std)
%MONTECARLO randomizes input variables along their provided standard
%deviations
%Inputs: Default parameters and standard deviation vector
%Outputs: Varied parameters

%Varies all parameters according to their standard deviation in the sigma
%array. Some varied by percent, others absolut magnitude.
F_T = F_T*(1-std(1)/100*randn(1, 1));
Mass_Initial = Mass_Initial*(1-std(2)/100*randn(1, 1));
Mass_Final = Mass_Final*(1-std(3)/100*randn(1, 1));
Launcher_Azimuth = Launcher_Azimuth+std(4)*randn(1, 1);
Launcher_Angle = Launcher_Angle+abs(std(5)*randn(1, 1));
C_b = C_b*(1-std(8)/100*randn(1, 1));
C_c = C_c*(1-std(9)/100*randn(1, 1));
F_Cant = F_Cant+std(10)*randn(1, 1);
F_T_a = F_T_a+std(11)*randn(1, 1);
F_T_o = F_T_o+std(12)*rand(1, 1);

end