%% plot proportion noise v. proportion convergence

figure % new figure
propCon = propCon_S0;
% specify colormap
P = [1, 2, 3, 6, 7, 10]; % noise vector
colormap(bone(length(P)));
customColor = bone(length(P)+2);
customColor = customColor(1:end-1,:); % remove white

for NP=1:length(P)
    plot(propCon{1,NP}, 'LineWidth', 2.5, 'Color', customColor(NP, :));
    hold on
end

ax = gca;
ax.FontSize = 25; 
ax.XTick = [1 2 3 4 5 6 7 8 9 10 11];
ax.XTickLabel = {'0' '.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9', '1'};
yticks([0 .5 1])
ylim([0, 1])
title("Network Capacity:S0", 'Fontsize', 30, 'fontname', 'calibri')
ylabel("Proportion convergence", 'Fontsize', 30,'fontname', 'calibri')
xlabel("Proportion noise", 'Fontsize', 30,'fontname', 'calibri')
lgd = legend('.02', '.04', '.06', '.12', '.14*', '.2'); % change this to alpha value
lgd.Title.String = 'alpha = P/N';

%% determine proportion of times that the network converged for different levels of noise

% for col=1:length(minCon(1,:))
%     propConverged(col) = 1 - (sum(isnan(minCon(:,col)))/length(minCon(:,col)));
% end