function EYt = computeEYt(x, ptm1, Yt, mu, sigma_dt)
%   ptm1 - discretized distribution of x at time t - 1 (point masses)
%   drift - a + b E^{\mu_{t-1}^X}[Y_{t-1}] + c Y_{t-1}: should be a vector
%   of length corresponding to the gridpoints of x

%   Using repmat representation from before

    integrand = ComputeYt(x, Yt, mu, sigma_dt);

    EYt = sum(ptm1 .* integrand);
end

