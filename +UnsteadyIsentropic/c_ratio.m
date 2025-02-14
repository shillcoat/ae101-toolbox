function ratio = c_ratio(du_c1,gamma)
    %C_RATIO Summary of this function goes here
    %   For a left-facing expansion wave, must set du_c1<0
    ratio = 1-((gamma-1)./2).*du_c1;
end

