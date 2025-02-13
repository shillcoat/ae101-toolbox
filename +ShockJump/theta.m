function t = theta(M1,gamma,beta)
    %THETA Deflection angle of oblique shock
    %   Returns oblique shock deflection angle
    M = abs(M1);  % Make sure Mach number is positive
    t = atan(cot(beta).*(M.^2.*(sin(beta)).^2-1)./...
        ((gamma+1)./2.*M.^2-(M.^2.*(sin(beta)).^2-1)));
end