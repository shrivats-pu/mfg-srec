function  yt ...
    =  StepBack_given_mu(x, pt, yt, ytp1,  dt, ...
                          iter, a, b, c, sigma, pi_k, gamma, K)
    % x - grid of banked SRECs
    % pt - mean field distribution across each agent type at time t (passed
    % through as squeezed, so it's only two dimensional (state and agent
    % space))
    % yt - estimate of adjoint of FBSDE at time t (passed through as
    % squeezed, so only two dimensional (state and agent space))
    % ytp1 - estimate of adjoing of FBSDE at time t+1 (passed through as
    % squeezed, so only two dimensional (state and agent space))
    % dt - increment of time-steps
    % iter - prespecified amount of internal iterations to use here 
    % a, b, c - constants (which may vary by agent type) that are used in
    % the forward equation of the FBSDE
    % sigma - volatility of agent types
    % pi_k - proportion of agent types
    % gamma - vector of trading cost parameters
    % K - number of agent types
                      
    x_r = repmat(x,1,length(x));
    
    % now iterate across iter
    for m = 1 : iter
        % calculation of E^{mu_t^k}[Y_t^k] (results in a number)
        EYt_all = sum( (pi_k ./ gamma) .* sum(  pt(:, :).* yt(:,:), 1) ) ;

        % iterate across agent types
        for k = 1:K
            % defining a_k, b_k, c_k for each agent type
            a_k = a(k);
            b_k = b(k);
            c_k = c(k);
            
            sigma_dt_k = sigma(k)*sqrt(dt); % variance of W

            ytp1_r = repmat(ytp1(:,k),1,length(x));

            mu_dt = x + (a_k + b_k* EYt_all + c_k * yt(:,k)  )* dt; % mean of X_{t+1}

            mu_dt_r = repmat(mu_dt',length(x),1); % repeating it for vectorization
            
            % computes Yt = E[Y_{t+1}|F_t] given the initial points of x 
            % this is the expectation with respect to the randomness of the
            % brownian motion (along paths of x across time, as opposed to along
            % the distribution of x for a frozen point in time (across agents)
            % does this in a vectorized manner for all grid-points of x
            % this is possible because of a locally linear assumption
            yt(:,k) = ComputeYt(x_r, ytp1_r, mu_dt_r, sigma_dt_k);
        end
    end
    
end

