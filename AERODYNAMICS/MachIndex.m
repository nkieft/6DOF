function [C_index] = MachIndex(M, Len)

M_index = round(M/0.01);

if(M_index > Len)
    M_index = Len-1;
    disp("WARNING: SIMULATION MAY HAVE DIVERGED");
end

C_index = M_index+1;

end