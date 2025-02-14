function ratio = p_ratio(du_c1,gamma)
    %P_RATIO Summary of this function goes here
    %   For a left-facing expansion wave, must set du_c1<0
    ratio = (1-((gamma-1)./2).*du_c1).^(2.*gamma./(gamma-1));
end

