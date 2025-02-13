function ratio = rho_ratio(M1,M2,gamma)
    %RHO_RATIO Density ratio relation for steady, isentropic flow
    %   Returns the steady isentropic density ratio rho2/rho1
    ratio = ((1+(gamma-1)./2.*M1.^2).^(1/(gamma-1)))./...
        ((1+(gamma-1)./2.*M2.^2).^(1/(gamma-1)));
end