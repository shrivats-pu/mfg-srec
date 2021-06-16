function S_finite =  FinitePlayerGame(x, p, Y, zeta, gamma, sigma, h,R, P, t, K, pi_k, total_players, mu_0, s_0)
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
    
    v = VideoWriter('movie.avi');
    open(v);
    for n = 1 : size(Y,2) 
        % plotting the mean field and empirical distribution 
        Plot_MFD_ED(bank(n,:), x, pi_k, p(:,n,:) );
        frame = getframe(gcf);
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
    close(v);

    % plotting the implied SREC price in the finite player game
    getframe()
    f1 = figure(21);
    ax = gca;
    ax.FontSize = 14;
    plot(t, S_sim)
    xlabel("Time", 'fontsize', 14)
    ylabel("SREC Price", 'fontsize', 14)
    title("Finite player game - implied SREC price", 'fontsize',14)

    names = {'Sub-population 1','Sub-population 2'};
    
    % plotting the controls for every agent across time
    getframe();
    for k = 1:1:K
       f2 = figure(101);
       ax = gca;
       ax.FontSize = 14;
       idx_agent_k = (agent_type == k);
       linS = {'-','-',':'};
       col = {[1, 0, 0, 0.5], [0, 0, 1, 0.5]};
       plot(t, g_sim(:, idx_agent_k), 'linestyle',linS{k}, 'color', col{k})
       title("Finite player game - optimal firm generation", 'fontsize', 14)
       xlabel("Time", 'fontsize', 14)
       ylabel("Generation Rate", 'fontsize', 14)
       %legend(names)
       hold on
       f3 = figure(102);
       ax = gca;
       ax.FontSize = 14;
       plot(t, Gamma_sim(:, idx_agent_k), 'linestyle',linS{k}, 'color', col{k})
       title("Finite player game - optimal firm trading", 'fontsize', 14)
       xlabel("Time", 'fontsize', 14)
       ylabel("Trading Rate", 'fontsize', 14)
       %legend(names)
       hold on
       
    end
    
    % plotting histograms of total banked SRECs, total generation, total
    % trading - note this will have to reworked for K > 2 
    f4 = figure(201);
    x0=10;
    y0=10;
    width=1200;
    height=width;
    set(gcf,'position',[x0,y0,width,height])
    for k = K:-1:1
        % initial / terminal SRECs
        subplot(2, 3, 1 + 3 * (k-1))
        ax = gca;
        ax.FontSize = 14;
        idx_agent_k = (agent_type == k);
        histogram(bank(end, idx_agent_k), x, 'edgecolor', 'None', 'facealpha', 1)
        hold on 
        histogram(bank(1, idx_agent_k), x, 'edgecolor', 'None', 'facealpha', 1)
        xlim([-0.1, 1.2])
        xline(R, ':r', 'LineWidth', 2.2)
        xlabel("Banked SRECs", 'fontsize', 14)
        ylabel("Frequency", 'fontsize', 14)
        title(sprintf('Initial / Terminal SRECs for Sub-population %d', k), 'fontsize', 14)
        
        % SREC generation (we sum until end-1 as there is no action taken
        % at t = T)
        subplot(2, 3, 2 + 3 * (k - 1))
        ax = gca;
        ax.FontSize = 14;
        histogram(sum(g_sim(1:end - 1, idx_agent_k), 1) * dt, x, 'edgecolor', 'None', 'facealpha', 1)
        xline(mean(sum(g_sim(1:end - 1, idx_agent_k), 1) * dt), ':r', 'LineWidth', 2.2)
        xlim([0.25 1.15]) % hard coded limit
        xlabel("Generated SRECs", 'fontsize', 14)
        ylabel("Frequency", 'fontsize', 14)
        title(sprintf('Generated SRECs for Sub-population %d', k), 'fontsize', 14)
        
        % SREC trading (we sum until end - 1 as there is no action taken at
        % t = T)
        subplot(2, 3, 3 + 3 * (k - 1))
        ax = gca;
        ax.FontSize = 14;
        histogram(sum(Gamma_sim(1:end-1, idx_agent_k), 1) * dt, x, 'edgecolor', 'None', 'facealpha', 1)
        xline(mean(sum(Gamma_sim(1:end-1, idx_agent_k), 1) * dt), ':r', 'LineWidth', 2.2)
        xlim([-0.25 0.25]) % hard coded limit
        xlabel("Traded SRECs", 'fontsize', 14)
        ylabel("Frequency", 'fontsize', 14)
        title(sprintf('Traded SRECs for Sub-population %d', k), 'fontsize', 14)
    end
    
    for k = 1:1:K
        % plotting correlation between total generation and total trading
        f5 = figure(301);
        ax = gca;
        ax.FontSize = 14;
        idx_agent_k = (agent_type == k);
        col = ['r', 'b'];
        scatter(sum(g_sim(1:end - 1, idx_agent_k), 1) * dt, sum(Gamma_sim(1:end-1, idx_agent_k), 1) * dt, col(k), '.')
        xlabel("Generated SRECs", 'fontsize', 14)
        ylabel("Traded SRECs", 'fontsize', 14)
        title("Generated SRECs vs Traded SRECs", 'fontsize', 14)
        %legend(names)
        hold on
    end

    % right now, this results in a perfectly linear relationship between
    % them - hypothesized that the slope of the line is related to the
    % ratio of zeta and gamma
    S_finite = S_sim;
    
    non_compliance_prob = zeros(K,1);
    for k = 1:K
        idx_agent_k = (agent_type == k);
        non_compliance_prob(k) = sum(bank(end, idx_agent_k) < R) / (pi_k(k)*total_players);
    end
    disp(non_compliance_prob);
%     save2pdf("finite_player_SREC.pdf", f1, 600)
%     save2pdf("finite_player_generation.pdf", f2, 600)
%     save2pdf("finite_player_trading.pdf", f3, 600)
%     save2pdf("finite_player_histograms.pdf", f4, 600)
%     save2pdf("gen_trade_relationship.pdf", f5, 600)
%     close all
%     save('finite_player_game.mat')
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
                
