% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 8, 2020

pattern = 2;
noise_iter = 9;
count = 1
figure
for i = 126:150
subplot(5,5,count)
plot(Convergence{1,pattern}{i, noise_iter}, 'LineWidth',1)
hold on
plot(m2Convergence{1,pattern}{i, noise_iter}, 'LineWidth',1)
count = count+1
end
