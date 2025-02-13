function ratio = p_ratio(M1,M2,gamma)
    %P_RATIO Pressure ratio relation for steady, isentropic flow
    %   Returns the steady isentropic pressure ratio p2/p1
    ratio = ((1+(gamma-1)./2.*M1.^2).^(gamma/(gamma-1)))./...
        ((1+(gamma-1)./2.*M2.^2).^(gamma/(gamma-1)));
end