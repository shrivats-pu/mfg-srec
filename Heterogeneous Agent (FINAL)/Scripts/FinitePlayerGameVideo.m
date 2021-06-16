function S_finite =  FinitePlayerGameVideo(x, p, Y, zeta, gamma, sigma, h,R, P, t, K, pi_k, total_players, mu_0, s_0)
    % x - vector of gridpoints 
    % p - theoretical mean field distribution of agents
    % Y - mean field solution of non-compliance probability for agents of
    % each type, at each value of x
    % zeta - vector of generation cost parameters (length K)
    % gamma - vector of trading cost parameters (length K)
    % sigma - vector of volatilities (length K)
    % h - vector of agents baseline generation rate (length K)
    % R - SREC requirement (scalar)
    % P - SACP (scalar)
    % t - vector of time-steps at which firms make decisions
    % K - number of agent classes
    % pi_k - long-run proportion of agents within each class
    % total_players - total number of agents in the system
    % mu_0 - mean of initial distribution
    % s_0 - standard deviation of initial distribution
    
    % This function takes in the above arguments and simulates a compliance
    % period forward for an SREC market with total_players agents - the
    % agents behave based on the simulated price (not the deterministic
    % mean field equilibrium SREC price). 
    
    
    % Potential future ideas to add to this: 
        % 1) Summary statistics and histograms / density plots for total g,
        % total Gamma, ending amount of banked SRECs
 
    
    %close all
    
    num_players = pi_k*total_players; % num_players defined
    dt = t(2) - t(1);
    y1 = ceil(total_players*0.025);
    y2 = ceil(total_players*0.05);
    y3 = ceil(total_players*0.075);

    % defining objects that will store data
    S_sim = NaN(size(Y,2),1);
    S_mfg = NaN(size(Y,2),1);
    g_sim = NaN(size(Y,2), total_players);
    Gamma_sim = NaN(size(Y,2), total_players);
    bank = NaN(size(Y,2), total_players);
    
    agent_type = repelem(1:K, num_players); % creating a list that tells us what type each agent is in
    % first N_k_1 agents of type k_1, next N_k_2 agents of type k_2
    
    % setting the first element of bank to be drawn from the initial
    % distribution for each type of agent
    Y_mean = zeros(1,K);
    for k = 1:K
        bank(1, agent_type == k) = mu_0(k) + s_0(k) * randn(sum(agent_type == k), 1);
    end
    f6 = figure(235);
    x0=10;
    y0=10;
    width=650;
    height=width;
    v = VideoWriter('1000_agents.avi');
    v.FrameRate = 10;
    open(v);
    for m = 1: 1%K+1
        for n = 1 : size(Y,2) 
            % plotting the mean field and empirical distribution
            %Plot_MFD_ED(bank(n,:), x, pi_k, p(:,n,:) );
            pp = p(:,n,:);
            % 1st subplot
            f6a = subplot(K+1, 1, m);
            
            h1 = histogram(bank(n,agent_type == 1), x, 'edgecolor','none');
            hold on;
            %f = zeros(length(x),1);
            f =  squeeze(pp(:,1,1));
            p1 = plot(x, f*length(bank(n,agent_type == 1)), 'linewidth',2  );
            ax = gca;
            ax.FontSize = 14;
            xlabel('SREC Inventory')
            ylabel("Density")
            ylim([0 y1])
            xlim([0 1.2])
            %title("Mean Field Distribution for Sub-population 1")
            hold off;
            % 2nd subplot
            f6b = subplot(K+1, 1, 2);
            h2 = histogram(bank(n,agent_type == 2), x, 'edgecolor','none');
            hold on;
            f = squeeze(pp(:,1,2));
            p2 = plot(x, f*length(bank(n,agent_type == 2)), 'linewidth',2  );
            ax = gca;
            ax.FontSize = 14;
            xlabel('SREC Inventory')
            ylabel("Density")
            ylim([0 y2])
            xlim([0 1.2])

            %title("Mean Field Distribution for Sub-population 2")
            hold off;
            % 3rd subplot
            f6c = subplot(K+1, 1, 3);
            
            h2 = histogram(bank(n,:), x, 'edgecolor','none');
            hold on;
            f = zeros(length(x),1);
            for k = 1 : size(pp,3)
                f = f + pi_k(k) * squeeze(pp(:,1,k));
            end
            p2 = plot(x, f*length(bank(n,:)), 'linewidth',2  );
            ax = gca;
            ax.FontSize = 14;
            xlabel('SREC Inventory')
            ylabel("Density")
            ylim([0 y3])
            xlim([0 1.2])

            %title("Mean Field Distribution for All Agents")
            hold off
            frame = getframe(f6);
            writeVideo(v,frame);
            for k = 1:K
                % updating Y_k for the banked values for each of the finite
                % number of firms
                x_k = bank(n, agent_type == k);
                Y_k = interp1(x, Y(:,n,k), x_k, 'linear', 'extrap');

                % if we're at the terminal condition, then Y is just the
                % indicator function
                if n == size(Y,2)                
                    Y_k = x_k < R;
                end
                % using the empirical mean instead of the expectation of Y 
                Y_mean(k) = mean(Y_k); 

            end

            % using the finite-player formula for equilibrium SREC price
            S_sim(n) = sum(P ./ gamma .* pi_k .* Y_mean) / sum(pi_k ./ gamma);

            % calculating the equilibrium SREC price using the mean field
            % formula
            numer = 0;
            for k = 1 : K
                numer = numer + sum( P /gamma(k) * pi_k(k) * sum(p(:,n,k) .* Y(:,n,k) ) );
            end
            S_mfg(n) = numer / sum(pi_k ./ gamma);                

            % the block of code below calculates the optimal behaviours for
            % each of the agents and updates their states
            if n < size(Y,2)
                for k = 1: K

                    idx_agent_k = (agent_type == k); % indicator to stratify behaviour by agent type

                    x_k = bank(n, idx_agent_k);
                    Y_k = interp1(x, Y(:,n,k), x_k, 'linear', 'extrap');

                    % calculation of the agents optimal behaviours per the
                    % formulas arising from the variational analysis
                    g_sim(n, idx_agent_k) = h(k) + P / zeta(k) * Y_k;
                    Gamma_sim(n, idx_agent_k) = 1 / gamma(k) * (P * Y_k - S_sim(n));

                    % updating the banked SRECSe of each firm
                    bank(n+1, idx_agent_k) = bank(n, idx_agent_k) + ...
                        (g_sim(n, idx_agent_k) + Gamma_sim(n, idx_agent_k)) * dt ...
                        + sigma(k) * randn(1, length(x_k)) * sqrt(dt);                

                end
            end
        end
    end
    close(v);
    S_finite = 0;
    %close all
end



% possibly useful debugging code if necessary
%                 g_str(n, idx_agent_k) = interp1(x, g(:, n,k), x_k,'linear', 'extrap');
%                 Gamma_str(n, agent_type == k)= interp1(x, Gamma(:,n,k),x_k,'linear', 'extrap');
% 
%                 g_diff(n, idx_agent_k) = g_str(n, idx_agent_k) - g_sim(n, idx_agent_k);
%                 Gamma_diff(n, idx_agent_k) = Gamma_str(n, idx_agent_k) - Gamma_sim(n, idx_agent_k);
% 
%                 bank(n+1, agent_type == k) = bank(n, idx_agent_k) + ...
%                     (g_str(n, idx_agent_k) + Gamma_str(n,idx_agent_k)) * dt ...
%                     + sigma(k) * randn(1, length(x_k))*sqrt(dt);
                
