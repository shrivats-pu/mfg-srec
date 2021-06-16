load('finite_player_game.mat')
load('base_scenario_SREC_MFG.mat')
f6 = figure(235);
x0=10;
y0=10;
width=600;
height=width;
pp = zeros(length(x), length(t));
for k = 1:K
    pp = pp + pi_k(k) .* p(:,:,k);
end
set(gcf,'position',[x0,y0,width,height])
Nplayers = 2000;
subset = [1, 14, 27, 40, 53];
col = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
   
for k = 1:K
    subplot(K+1, 1, k)
    idx_agent_k = (agent_type == k);
    for a = 1:length(subset)
        plot(x, p(:, subset(a), k)*Nplayers*pi_k(k), 'color', col{a});
        hold on
        histogram(bank(subset(a), idx_agent_k),x, 'EdgeColor', 'none', 'EdgeAlpha', 0, 'FaceAlpha', 0.2, 'FaceColor',col{a});
    end
    ax = gca;
    ax.FontSize = 14;
    xlabel('SREC Inventory')
    xlim([0 max(x)])
    ylabel("Density")
    title(strcat("Mean Field Distribution for Sub-population ", int2str(k)))
end

subplot(K+1, 1, K+1)
plot(x, pp(:, subset)*Nplayers);
hold on
for b = 1:length(subset)
    histogram(bank(subset(b), :), x, 'EdgeColor', 'none', 'EdgeAlpha', 0, 'FaceAlpha', 0.2, 'FaceColor',col{b})
end
ax = gca;
ax.FontSize = 14;
xlabel('SREC Inventory');
ylabel("Density")
title("Mean Field Distribution for All Agents")
xlim([0 max(x)])
save2pdf("finite_player_dist.pdf", f6, 600)
