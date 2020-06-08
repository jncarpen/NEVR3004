% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 8, 2020

Convergence = m2Convergence;
proportionNoise = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
c = lines(20); % set colormap
for pattern = 2
    % Recall: P = [1, 2, 3, 6, 7, 10];
    figure % intialize new figure
    set(gcf,'color','w'); % set background color to white
    for noise_iter = 2:11 % loop through all possible noise values
        
        % create a subplot for each noise category
        subplot(2,5,noise_iter-1)

        for i=1:50:1000
            plot(Convergence{1,pattern}{i,noise_iter}, 'LineWidth', 1)
            ax = gca;
            ax.FontSize = 17; 
            yticks([-1 0 1])
            xticks([100, 300, 500])
            xlim([1,500]);
            ylim([-1, 1])
            hold on
        end
        
        caption = sprintf('%d %% noise', proportionNoise(noise_iter)*100);
        title(caption, 'FontSize', 18)
        
    end
end
