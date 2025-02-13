function mach = M2(M1,gamma,beta,t)
    %M2 Mach number downstream of shock jump condition
    %   Returns the post-shock Mach number M2
    M = abs(M1);  % Make sure Mach number is positive
    mach = sqrt((1+(gamma-1)./(gamma+1).*((M.*sin(beta)).^2-1))./...
        (1+(2.*gamma)./(gamma+1).*((M.*sin(beta)).^2-1))./...
        (sin(beta-t).^2));
end