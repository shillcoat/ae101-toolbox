function ratio = h_ratio(M1,gamma,beta)
    %H_RATIO Enthalpy ratio relation for shock jump condition
    %   Returns the shock jump condition enthalpy ratio h2/h1=T2/T1
    M = abs(M1);  % Make sure Mach number is positive
    ratio = 1+2.*(gamma-1)./(gamma+1).^2.*...
        ((M.*sin(beta)).^2-1).*(1+gamma.*(M.*sin(beta)).^2)./...
        ((M.*sin(beta)).^2);
end