% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 5, 2020

figure
for i=1:10:1000
plot(ConvergenceData{1,3}{i,3}, 'LineWidth', 1.25)
hold on
end