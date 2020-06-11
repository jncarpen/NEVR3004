% make sample trajectory plots for 2 patterns 
% Last modified: June 8, 2020

Convergence = Convergence_NS2;
m2Convergence =m2Convergence_NS2;
pattern = 2;
% Recall: P = [1, 2, 3, 6, 7, 10];
noise_iter = 10;
count = 1;
figure
for i = 1
% subplot(5,5,count)
% plot overlap of state with cued pattern
plot(Convergence{1,pattern}{i, noise_iter}, 'LineWidth',1) 
hold on
% plot overlap of state with distractor pattern
plot(m2Convergence{1,pattern}{i, noise_iter}, 'LineWidth',1, 'color', 'black')
ax = gca;
ax.FontSize = 20; 
yticks([-1 0 1])
xticks([100, 300, 500])
xlim([1,500]);
ylim([-1, 1])
xlabel("noise", 'Fontsize', 30,'fontname', 'calibri')
ylabel("convergence", 'Fontsize', 30,'fontname', 'calibri')
title(sprintf('random(-):iter%d', i))
lgd = legend('cued pattern', 'distractor pattern');
count = count+1;
end
