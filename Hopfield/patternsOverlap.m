% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 8, 2020

Convergence = Convergence_NS2;
m2Convergence =m2Convergence_NS2;
pattern = 5;
% Recall: P = [1, 2, 3, 6, 7, 10];
noise_iter = 10;
count = 1;
figure
for i = 51:75
subplot(5,5,count)
plot(Convergence{1,pattern}{i, noise_iter}, 'LineWidth',2)
hold on
plot(m2Convergence{1,pattern}{i, noise_iter}, 'LineWidth',2, 'color', 'black')
% ax = gca;
% ax.FontSize = 20; 
% yticks([-1 0 1])
% xticks([100, 300, 500])
% xlim([1,500]);
% ylim([-1, 1])
% xlabel("proportion noise", 'Fontsize', 30,'fontname', 'calibri')
% ylabel("proportion convergence", 'Fontsize', 30,'fontname', 'calibri')
title(sprintf('NS0:iter%d', i))
% lgd = legend('Pattern 1', 'Pattern 2');
count = count+1;
end
