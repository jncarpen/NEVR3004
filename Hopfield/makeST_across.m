% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 8, 2020

Convergence = Convergence;
proportionNoise = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
c = lines(20); % set colormap
for pattern = 2
    % Recall: P = [1, 2, 3, 6, 7, 10];
    figure % intialize new figure
    for noise_iter = 3 % loop through all possible noise values
        
        % create a subplot for each noise category
%         subplot(3,4,noise_iter)
%         title(sprintf('%d percent noise',proportionNoise(noise_iter)*100))
%         set(gcf,'color','w'); % set background color to white

        for i=1:20:1000
            plot(Convergence{1,pattern}{i,noise_iter}, 'LineWidth', 1)
            ax = gca;
            ax.FontSize = 30; 
            yticks([-1 0 1])
            xticks([100, 300, 500])
            xlim([1,500]);
            ylim([-1, 1])
            hold on
        end
    end
end
