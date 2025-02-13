function ratio = rho_ratio(M1,gamma,beta)
    %RHO_RATIO Density ratio relation for shock jump condition
    %   Returns the shock jump condition density ratio rho2/rho1
    M = abs(M1);  % Make sure Mach number is positive
    ratio = (M.*sin(beta)).^2./(1+(gamma-1)./(gamma+1).*...
        ((M.*sin(beta)).^2-1));
end