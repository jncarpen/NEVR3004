%% Asynchronous non-sequential update 

N = 50; % size of network
[V, W] = patternWeight(N); % generate random pattern & weights matrix

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
figure() % initialize figure

for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
    S = V; % initialize network
    x = addNoise(S, N, proportionNoise(noise_iter)); % create vector w/ varied amounts of noise
    m(1) = (x * V')/N; % calculate overlap between state & pattern
    
    overlap = []; % initialize an overlap matrix
       
        for simLength = 1:300 % total # of iterations *
            
            neuron = randi([1 50]); % random number from 1-50
            h(neuron) = (W(neuron,:) * x')/N; % compute input potential of neuron_i
            x(neuron) = sign(h(neuron)); % update the state of the network
            m(neuron) = (x*V')/N; % similarity between current state & pattern
            overlap = [overlap, m(neuron)]; %store in overlap matrix
            
        end
        
        subplot(5,2,noise_iter)
        plot(overlap)
        title(sprintf('%dth percent noise',proportionNoise(noise_iter)))
end