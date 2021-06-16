function Plot_MFD_ED(x_emp, x_grid, pi, p)
    
    fig = figure(1);
    clf(fig);
    % plot empirical histogram distribution
    histogram(x_emp, x_grid, 'edgecolor','none');
    
    % figure labels
    xlabel("SREC Inventory", 'fontsize', 14)
    ylabel("Frequency", 'fontsize', 14)
    title("Empirical and Theoretical Distribution of Agents", 'fontsize', 14)
    
    hold on;
    
    % theoretical distribution (from solution to MFG)
    f = zeros(length(x_grid),1);
    for k = 1 : size(p,3)
        f = f + pi(k) * squeeze(p(:,1,k));
    end
    % scaled by the length of the empirical x_grid
    set(gca,'fontsize',16);
    xlabel('$x$','fontsize',20,'interpreter','latex');
    ylabel('$F(x)$','fontsize',20,'interpreter','latex');
    ylim([0 14]); %hard-coded for now
    plot(x_grid, f*length(x_emp), 'linewidth',2  );
    

end