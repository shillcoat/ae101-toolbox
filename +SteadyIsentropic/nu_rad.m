function n = nu_rad(M,gamma)
    %NU Prandtl-Meyer function in radians for steady, isentropic flow
    %   Returns the Prandtl-Meyer function for steady, isentropic flow.
    %   Note that units are radians.
    n = sqrt((gamma+1)./(gamma-1)).*atan(sqrt((gamma-1)./...
        (gamma+1).*(M.^2-1)))-atan(sqrt(M.^2-1));
end

