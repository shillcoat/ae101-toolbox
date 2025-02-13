function ratio = du_c1(M1,gamma)
    %DU_C1 Ratio of change in velocity to upstream speed of sound
    %   Returns the change in velocity across the shock, normalized by the
    %   upstream speed of sound du/c1
    M = abs(M1);  % Make sure Mach number is positive
    ratio = 2./(gamma+1).*(M.^2-1)./M;
end