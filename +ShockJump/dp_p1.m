function ratio = dp_p1(du_c1,gamma)
    %DP_P1 Ratio of change in pressure to upstream pressure
    %   Returns the change in pressure across the shock, normalized by the
    %   upstream pressure dp/p1
    ratio = gamma.*du_c1.*((gamma+1)./4.*du_c1+...
        sqrt(1+((gamma+1)./4).^2.*du_c1.^2));
end
