function [bank, cash, gen, trade] = SimPath(t, x, x0, g, Gamma, S, ...
                zeta, gamma, sigma, h, P, R, K, Y, NSims)
    % plots a sample path for a representative agent of each sub-population
    % through the compliance period
    % NSims should be 1 if you want the plots. If you want summary
    % statistics, just comment out the plotting. This was not used in the
    % paper, hence the shoddy commenting / hacky setup.
    dt = t(2)-t(1);
    bank = NaN(K, length(t), NSims);
    bank_planned = NaN(K, length(t), NSims);
    cash = NaN(K, length(t), NSims);
    gen = NaN(K, length(t), NSims);
    trade = NaN(K, length(t), NSims);
    rndmnss = NaN(K, length(t), NSims);
    non_compliance_prob = NaN(K, length(t), NSims);
    for m = 1:NSims
        bank(:,1, m) = x0;
        bank_planned(:,1, m)=x0;
        cash(:,1, m) = 0;

        for n = 1 : length(t)-1
            %W_t = randn(1,1); % same randomness to both firms
            for k = 1:K
                non_compliance_prob(k, n, m) = interp1(x,Y(:,n,k),bank(k,n, m), 'linear', 'extrap');

                gen(k,n, m) = interp1(x, g(:,n,k), bank(k,n, m), 'linear', 'extrap');
                trade(k,n, m) = interp1(x, Gamma(:,n,k), bank(k,n, m), 'linear', 'extrap');
                W_t = randn(1,1); % different randomness to both firms
                bank(k,n+1, m) = bank(k,n, m) ...
                    + (gen(k,n, m) + trade(k,n, m))*dt + sigma(k).*W_t*sqrt(dt);

                bank_planned(k,n+1,m) = bank(k,n,m) + (gen(k,n,m) + trade(k,n,m))*dt;        

                rndmnss(k, n+1,m) =  bank(k, n+1,m) - bank_planned(k, n+1,m);

                cash(k,n+1,m) = cash(k,n,m) ...
                    - (0.5*zeta(k)*(gen(k,n,m)-h(k))^2 + 0.5*gamma(k)*trade(k,n,m)^2 ...
                        +  S(n)*trade(k,n,m))*dt;
            end
        end

        non_compliance_prob(:,end,m) = bank(:,end,m) < R;
        cash(:,end,m) = cash(:,end,m) - P*max(R - bank(:,end,m), 0);

        fig = figure(30);
        clf(fig);
        x0=10;
        y0=10;
        width=700;
        height=width;
        set(gcf,'position',[x0,y0,width,height])
        subplot(6,1,1);
        col = {'r', 'b'};
        for k=1:K
            plot(t,bank(k,:), 'color', col{k});
            hold on
        end
        ax = gca;
        ax.FontSize = 17;
        ylabel('SRECs');
        title("Sample Path for Representative Firms in Mean Field Game")
        %plot(t, bank_planned(:,:), '--');
        %legend('Firm in Sub-population 1', 'Firm in Sub-population 2') 
        
        subplot(6,1,2);
        for k = 1:K
            plot(t,gen(k,:), 'color', col{k});
            hold on
        end
        ax = gca;
        ax.FontSize = 17;
        ylabel('Gen Rate');
        
        subplot(6,1,3);
        for k = 1:K
            plot(t,trade(k,:), 'color', col{k}); 
            hold on
        end
        ax = gca;
        ax.FontSize = 17;    
        ylabel('Trade Rate');
        
        subplot(6,1,4);
        for k = 1:K
            plot(t,cash(k,:), 'color', col{k});
            hold on
        end
        ax = gca;
        ax.FontSize = 17; 
        ylabel('Cum. Profit');
        
        subplot(6,1,5);
        for k = 1:K
            plot(t,cumsum(rndmnss(k,:), 'omitnan'), 'color', col{k});
            hold on
        end
        ax = gca;
        ax.FontSize = 17; 
        ylabel('Cum. Noise');
        
        subplot(6,1,6);
        for k = 1:K
            plot(t,non_compliance_prob(k,:), 'color', col{k});
            hold on
        end
        ax = gca;
        ax.FontSize = 17; 
        ylabel('$P(X_T < R)$','Interpreter','latex')
        xlabel('Time')
    %     
        save2pdf("path.pdf", fig, 600)
    end
    save('simulated_path_base.mat')
    terminal_SRECs_mean = NaN(K,1);
    terminal_SRECs_std = NaN(K,1);
    terminal_SREC_Q1 = NaN(K,1);
    terminal_SREC_Q3 = NaN(K,1);
    terminal_SREC_skew = NaN(K,1);
    terminal_SREC_kurt = NaN(K,1);
    
    planned_gen_SRECs_mean = NaN(K,1);
    planned_gen_SRECs_std = NaN(K,1);
    planned_gen_SREC_Q1 = NaN(K,1);
    planned_gen_SREC_Q3 = NaN(K,1);
    planned_gen_SREC_skew = NaN(K,1);
    planned_gen_SREC_kurt = NaN(K,1);
    
    trade_SRECs_mean = NaN(K,1);
    trade_SRECs_std = NaN(K,1);
    trade_SREC_Q1 = NaN(K,1);
    trade_SREC_Q3 = NaN(K,1);
    trade_SREC_skew = NaN(K,1);
    trade_SREC_kurt = NaN(K,1);
    
    profit_SRECs_mean = NaN(K,1);
    profit_SRECs_std = NaN(K,1);
    profit_SREC_Q1 = NaN(K,1);
    profit_SREC_Q3 = NaN(K,1);
    profit_SREC_skew = NaN(K,1);
    profit_SREC_kurt = NaN(K,1);
    
    non_compliance = NaN(K,1);
    % CALCULATE SUMMARY STATISTICS
    for k = 1:K
        terminal_SRECs_mean(k) = mean(bank(k, end, :));
        terminal_SRECs_std(k) = std(bank(k, end, :));
        terminal_SREC_Q1(k) = quantile(bank(k, end,:), 0.25);
        terminal_SREC_Q3(k) = quantile(bank(k, end,:), 0.75);
        terminal_SREC_skew(k) = skewness(bank(k, end, :));
        terminal_SREC_kurt(k) = kurtosis(bank(k, end, :));
        
        planned_gen_SRECs_mean(k) = mean(sum(gen(k, 1:end-1, :), 2)*dt);
        planned_gen_SRECs_std(k) = std(sum(gen(k, 1:end-1, :), 2)*dt);
        planned_gen_SREC_Q1(k) = quantile(sum(gen(k, 1:end-1, :), 2)*dt, 0.25);
        planned_gen_SREC_Q3(k) = quantile(sum(gen(k, 1:end-1, :), 2)*dt, 0.75);
        planned_gen_SREC_skew(k) = skewness(sum(gen(k, 1:end-1, :)*dt, 2));
        planned_gen_SREC_kurt(k) = kurtosis(sum(gen(k, 1:end-1, :)*dt, 2));

        trade_SRECs_mean(k) = mean(sum(trade(k, 1:end-1, :),2)*dt);
        trade_SRECs_std(k) = std(sum(trade(k, 1:end-1, :),2)*dt);
        trade_SREC_Q1(k) = quantile(sum(trade(k, 1:end-1, :), 2)*dt, 0.25);
        trade_SREC_Q3(k)= quantile(sum(trade(k, 1:end-1, :),2)*dt, 0.75);
        trade_SREC_skew(k) = skewness(sum(trade(k, 1:end-1, :)*dt,2));
        trade_SREC_kurt(k) = kurtosis(sum(trade(k, 1:end-1, :)*dt,2));

        profit_SRECs_mean(k) = mean(cash(k, end, :));
        profit_SRECs_std(k) = std(cash(k, end, :));
        profit_SREC_Q1(k) = quantile(cash(k, end, :), 0.25);
        profit_SREC_Q3(k) = quantile(cash(k, end, :), 0.75);
        profit_SREC_skew(k) = skewness(cash(k, end, :));
        profit_SREC_kurt(k) = kurtosis(cash(k, end, :));
        
        non_compliance(k) = sum(bank(k, end, :) < R) / NSims;

    end
    disp("TERMINAL SRECS")
    disp(terminal_SRECs_mean);
    disp(terminal_SRECs_std);
    disp(terminal_SREC_Q1);
    disp(terminal_SREC_Q3);
    disp(terminal_SREC_skew);
    disp(terminal_SREC_kurt);
    disp("GENERATED SRECS")
    disp(planned_gen_SRECs_mean);
    disp(planned_gen_SRECs_std);
    disp(planned_gen_SREC_Q1);
    disp(planned_gen_SREC_Q3);
    disp(planned_gen_SREC_skew);
    disp(planned_gen_SREC_kurt);
    disp("TRADED SRECS")
    disp(trade_SRECs_mean);
    disp(trade_SRECs_std);
    disp(trade_SREC_Q1);
    disp(trade_SREC_Q3);
    disp(trade_SREC_skew);
    disp(trade_SREC_kurt);
    disp("PROFIT")
    disp(profit_SRECs_mean);
    disp(profit_SRECs_std);
    disp(profit_SREC_Q1);
    disp(profit_SREC_Q3);
    disp(profit_SREC_skew);
    disp(profit_SREC_kurt);
    disp("NON-COMPLIANCE")
    disp(non_compliance)
    for k = K:-1:1
       figure(101+(k-1))
       linS = {'-','-',':'};
       col = {[0, 0, 1, 0.5], [1, 0, 0, 0.5]};
       for m = 1:1
            plot(t, bank(k, :,m), 'linestyle',linS{k}, 'color', col{k})
            hold on
       end
       
    end

    