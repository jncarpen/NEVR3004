% Correlation and confusion matrices

figure
% visualize patterns 
subplot(2,2,1)
imagesc(allPatterns);
xticks(0:10:50);
yticks(1:10);
title("pattern correlation",'Fontsize', 12,'fontname', 'calibri')
ylabel("pattern", 'Fontsize', 10, 'fontname', 'calibri')
xlabel("node",'Fontsize', 10,'fontname', 'calibri')
colorbar

% visualize correlation between patterns
subplot(2,2,2)
allPatterns(1,:); % pattern 1