% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 8, 2020

% Convergence = Convergence;
proportionNoise = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
for pattern = 1
    % Recall: P = [1, 2, 3, 6, 7, 10];
    figure % intialize new figure
    set(gcf,'color','w'); % set background color to white
    
    count = 1;
    for noise_iter = 2:2:11 % loop through all possible noise values
        % create a subplot for each noise category
        subplot(3,5,count)

        for i=1:30:1000
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
        count = count+1;
    end
end
