function ratio = c_ratio(M1,gamma,beta)
    %C_RATIO Speed of sound ratio relation for shock jump condition
    %   Returns the shock jump condition speed of sound ratio c2/c1
    M = abs(M1);  % Make sure Mach number is positive
    ratio = sqrt(1+2.*(gamma-1)./(gamma+1).^2.*...
        ((M.*sin(beta)).^2-1).*(1+gamma.*(M.*sin(beta)).^2)./...
        ((M.*sin(beta)).^2));
end