function p = UpdateMu(x, p0, Y, sigma, dt, a, b, c, gamma, pi_k)
    % x - the grid of points of banked SRECs
    % p0 - the initial distribution of banked SRECs across every agent
    % class
    % Y - the estimate of the adjoint process of the FBSDE (can be thought
    % of as the non-compliance probability)
    % sigma - volatilities for each type of agent
    % dt - length of timestep
    % a, b, c - constants (which may vary across agent types) that appear
    % in the forward equation of the FBSDE
    % gamma - vector of the trading cost parameters across agent types
    % pi_k - the proportion of agents of each type
    
    % generate the distribution X for each t given that
    % dX = ( a + b E[Y_t] + c Y_t ) dt + sigma dW_t
    % and X_0 ~ p0
    
    K = length(pi_k); % defining K to be the number of types of agents
    
    p = zeros(length(x), size(Y,2), K); % discretize distribution across a matrix for each x and every time-step

    p(:,1, :) = p0; 
    
    for n = 1 : size(Y,2)-1 % iterate across time-steps
        % calculate expected value of Y_t with respect to the distribution of the agents states
        EYt_all = sum( (pi_k ./ gamma) .* sum(  squeeze(p(:,n,:).* Y(:,n,:)), 1) ) ;
        % results in a number

        for k = 1:K
            % defining useful items for the FBSDE
            a_k = a(k);
            b_k = b(k);
            c_k = c(k);
            sigma_k = sigma(k);
            pp = zeros(length(x),1);

            parfor i = 1 : length(x)

                drift =  (a_k + b_k * EYt_all + c_k*Y(i,n,k))*dt; % updating drift for mean of X_{t+1}
                vol = sigma_k*sqrt(dt); % updating variance of X_{t+1}
                
                % because X_{t+1} | X_t is normal, we can use
                % GenDiscreteDist
                
                % essentially, this looks at every point in the grid of x
                % and determines the probability that you end up at said
                % point, conditional on the position of the agent at the
                % previous time-step. Summing this across every possible
                % starting point yields the probability that an agent falls
                % on a particular point at time t+1. The forward iteration
                % through time ensures we are always able to calculate this
                pp =  pp+ p(i,n,k) * GenDiscreteDist(x, x(i) + drift, vol);
                % the distribution is a mixture model 
            end

            p(:,n+1,k) = pp;% update the n+1th distribution 
        end
    end
    
    
end