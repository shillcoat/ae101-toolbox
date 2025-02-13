function ratio = h_ratio(M1,M2,gamma)
    %H_RATIO Enthalpy ratio relation for steady, isentropic flow
    %   Returns the steady isentropic enthalpy ratio h2/h1=T2/T1
    ratio = (1+(gamma-1)./2.*M1.^2)./(1+(gamma-1)./2.*M2.^2);
end