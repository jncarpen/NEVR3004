% make sample trajectory plots
% update this later
% Jo Carpenter
% Last modified: June 5, 2020

proportionNoise = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1];
for pattern = 1:6
    figure
    for noise_iter = 1:11 % loop through all possible noise values
        subplot(3,4,noise_iter)
        title(sprintf('%d percent noise',proportionNoise(noise_iter)*100))

        for i=1:50:1000
            plot(ConvergenceData{1,pattern}{i,noise_iter}, 'LineWidth', 1.25)
            xlim([1,500]);
            hold on
        end
    end
    %suptitle(sprintf('pattern=%u', pattern));
end