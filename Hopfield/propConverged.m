%% plot proportion noise v. proportion convergence

figure % new figure

% specify colormap
colormap(jet(P));
customColor = jet(P);

for NP=1:P
    plot(propCon2{1,NP}, 'LineWidth', 2, 'Color', customColor(NP, :));
    hold on
end

title("Network Capacity")
ylabel("Proportion convergence")
xlabel("Proportion noise")
lgd = legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'); % change this to alpha value
lgd.Title.String = 'Number of patterns';

%% determine proportion of times that the network converged for different levels of noise

% for col=1:length(minCon(1,:))
%     propConverged(col) = 1 - (sum(isnan(minCon(:,col)))/length(minCon(:,col)));
% end