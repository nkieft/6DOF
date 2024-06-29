function [C_b, C_c, Cp] = AeroImport(Max_Mach)
%AEROIMPORT Reads RASAero II "Run Test" Files to determine aerodynamic coefficients
%
%Inputs : RASAero II test files in directory
%
%Outputs: C_b, C_c, Cp
%C_b: Burn phase axial and normal force coefficients
%C_c: Coast phase axial and normal force coefficients
%Cp: Center of pressure (distance to tip) [m]
Max_Mach_index = round(Max_Mach*600);

%Loop through all files
for i = 1:12

%Parsing shenanigans
filename = strcat("CDDataFile", string(i-1), ".txt");
a = readmatrix(filename, "Delimiter", ["   ", "  "] );
a(isnan(a(:, 1)), :) = [];
a = a(1:Max_Mach_index, :);

%Extract data of interest from files
Cp(:, i) = a(4:6:end, 4);
CN_c(:, i) = a(5:6:end,5);
CN_b(:, i) = a(6:6:end, 5);
CA_c(:, i) = a(5:6:end, 6);
CA_b(:, i) = a(6:6:end, 6);

end

%Format data according to output
C_b(:, :, 1) = CN_b;
C_b(:, :, 2) = CA_b;
C_c(:, :, 1) = CN_c;
C_c(:, :, 2) = CA_c;

end