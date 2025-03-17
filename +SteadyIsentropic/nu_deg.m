function n = nu_deg(M,gamma)
    %NU Prandtl-Meyer function in degrees for steady, isentropic flow
    %   Returns the Prandtl-Meyer function for steady, isentropic flow.
    %   Note that units are degrees.
    n = sqrt((gamma+1)./(gamma-1)).*atand(sqrt((gamma-1)./...
        (gamma+1).*(M.^2-1)))-atand(sqrt(M.^2-1));
end

