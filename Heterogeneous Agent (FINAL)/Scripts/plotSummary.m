function plots = plotSummary(x, p, g, Gamma, S, h, t, K, pi_k)
% PLOT CONTROL SUMMARY_____________________________________________________
    y_max_g = max(g(:,:,:), [], 'all')*1.2;
    y_min_g = min(h)*0.8;

    y_max_Gamma = max(Gamma(:,:,:), [], 'all')*1.2;
    y_min_Gamma = min(Gamma(:,:,:), [], 'all')*1.2;
    f1 = figure(11);
    x0=10;
    y0=10;
    width=700;
    height=width;
    set(gcf,'position',[x0,y0,width,height])
    subset = [1, 14, 27, 40, 53];
    for k=1:K
        subplot(2,2,1 + (k-1)*2);
        plot(x,g(:,subset,k));
        ax = gca;
        ax.FontSize = 14;
        xlim([min(x) max(x)]);
        ylim([y_min_g y_max_g])
        %legend('$t = 0$', '$t = 0.25$' ,'$t = 0.5$', '$t = 0.75$', '$t = 1$','Interpreter','latex') % R2018b and later
        xlabel('Banked SRECs');
        ylabel("Generation Rate (SRECs / year)")
        title(strcat("Generation Rate for Sub-population ", int2str(k)))

        subplot(2,2,2 + (k-1)*2);
        plot(x,Gamma(:,subset,k));
        ax = gca;
        ax.FontSize = 14;
        xlim([min(x) max(x)]);
        ylim([y_min_Gamma y_max_Gamma])
        xlabel('Banked SRECs');
        ylabel("Trading Rate (SRECs / year)")
        title(strcat("Trading Rate for Sub-population ", int2str(k)))
        %legend('$t = 0$', '$t = 0.25$' ,'$t = 0.5$', '$t = 0.75$', '$t = 1$','Interpreter','latex') % R2018b and later

    end

    % PLOT EQUILIBRIUM PRICE___________________________________________________
    f2 = figure(12);
    x0=10;
    y0=10;
    width=600;
    height=width;
    set(gcf,'position',[x0,y0,width,height])
    plot(t, S)
    ax = gca;
    ax.FontSize = 14;
    xlabel("Time")
    ylabel("SREC Price")
    title("Equilibrium SREC Price", "FontSize", 14)
    disp(S)
    ylim([mean(S)*0.8 mean(S)*1.2])

    %PLOT MF DISTRIBUTION THROUGH TIME_________________________________________
    % presentation only
    %f3 = figure(13);
    pp = zeros(length(x), length(t));
    for k = 1:K
        pp = pp + pi_k(k) .* p(:,:,k);
    end
%     x0=10;
%     y0=10;
%     width=1200;
%     height=width*0.5;
%     set(gcf,'position',[x0,y0,width,height])
%     for l = 1:length(t)
%         subplot(1, 3, 1)
%         plot(x, p(:,l,1))
%         ax = gca;
%         ax.FontSize = 14;
%         xlim([min(x) max(x)])
%         ylim([min(p(:,:,1), [], 'all') max(p(:,:,1), [], 'all')*1.1])
%         ylabel("Density")
%         title("Mean Field Distribution for Sub-population 1")
%         subplot(1, 3, 2)
%         plot(x, p(:,l,2))
%         ax = gca;
%         ax.FontSize = 14;
%         xlim([min(x) max(x)])
%         ylim([min(p(:,:,2), [], 'all') max(p(:,:,2), [], 'all')*1.1])
%         xlabel("Banked SRECs")
%         title("Mean Field Distribution for Sub-population 2")
%         subplot(1, 3, 3)
%         plot(x, pp(:,l))
%         ax = gca;
%         ax.FontSize = 14;
%         xlim([min(x) max(x)])
%         ylim([min(pp(:,:), [], 'all') max(pp(:,:), [], 'all')*1.1])
%         title("Mean Field Distribution of Population")
%         
%     end
    
    f4 = figure(235);
    x0=10;
    y0=10;
    width=600;
    height=width;
    set(gcf,'position',[x0,y0,width,height])
    for k = 1:K
        subplot(K+1, 1, k)
        plot(x, p(:, subset, k));
        ax = gca;
        ax.FontSize = 14;
        xlabel('SREC Inventory')
        xlim([0 max(x)])
        ylabel("Density")
        title(strcat("Mean Field Distribution for Sub-population ", int2str(k)))
    end
    
    subplot(K+1, 1, K+1)
    plot(x, pp(:, subset));
    ax = gca;
    ax.FontSize = 14;
    xlabel('SREC Inventory');
    ylabel("Density")
    title("Mean Field Distribution for All Agents")
    xlim([0 max(x)])

    f = figure(234);
    for k = 1:K
        subplot(K+1,1,k)
        surf(t, x, p(:,:,k), 'EdgeColor', 'None')
        ax = gca;
        ax.FontSize = 14;
        xlabel("Time")
        ylabel("SREC Inventory")
        zlabel("Density")
        title(strcat("Mean field distribution for sub-population ", int2str(k)))
    end

    subplot(K+1,1,K+1)
    x0=10;
    y0=10;
    width=600;
    height=width;
    set(gcf,'position',[x0,y0,width,height])
    surf(t, x, pp, 'EdgeColor', 'None')
    ax = gca;
    ax.FontSize = 14;
    xlabel("Time")
    ylabel("SREC Inventory")
    zlabel("Density")
    title("Mean field distribution for all agents ")

    plots = [f1, f2, f4, f];
    %save2pdf("optimal_behaviour_base.pdf", f1, 600)
    %save2pdf("SREC_inf_pop.pdf", f2, 600)
    %save2pdf("mean_fields.pdf", f, 600)
    %save2pdf("mean_fields_2d.pdf", f4, 600)

end

