function ratio = rho_ratio(du_c1,gamma)
    %RHO_RATIO Summary of this function goes here
    %   For a left-facing expansion wave, must set du_c1<0
    ratio = (1+((gamma-1)./2).*du_c1).^(2./(gamma-1));
end

