function C = AeroData(C_in, AoA)
%AeroData Interpolates a vector of aero coefficients at a given mach for a
%specific angle of attack

    %High and low indices
    AoA_low = floor(AoA);
    AoA_high = ceil(AoA);
    index_low = AoA_low+1;
    index_high = AoA_high+1;

    %For high AoAs
    if(AoA_low > 10)
        AoA_low = 10;
        index_low = 11;
    end

    %For high AoAs
    if(AoA_high > 10)
        AoA_high = 45;
        index_high = 12;
    end

    %Interpolation
    y1 = C_in(index_low);
    y2 = C_in(index_high);
    x1 =AoA_low;
    x2 = AoA_high;

    C = y1+(AoA-x1)*(y2-y1)/(x2-x1);

    if(isnan(C))
        C = y1;
    end

    
end