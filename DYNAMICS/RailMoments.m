function M_R = RailMoments(onRail, M_F)

if(onRail)
    M_R = -M_F;
else
    M_R = zeros(1, 3);
end

end