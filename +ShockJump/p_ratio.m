function ratio = p_ratio(M1,gamma,beta)
    %P_RATIO Pressure ratio relation for shock jump condition
    %   Returns the shock jump condition pressure ratio p2/p1
    M = abs(M1);  % Make sure Mach number is positive
    ratio = 1+2.*gamma/(gamma+1).*((M.*sin(beta)).^2-1);
end