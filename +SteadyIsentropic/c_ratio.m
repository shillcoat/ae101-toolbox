function ratio = c_ratio(M1,M2,gamma)
    %C_RATIO Speed of sound ratio relation for steady, isentropic flow
    %   Returns the steady isentropic speed of sound ratio c2/c1
    ratio = ((1+(gamma-1)./2.*M1.^2).^(1/2))./...
        ((1+(gamma-1)./2.*M2.^2).^(1/2));
end