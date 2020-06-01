%% Task 2: Retrieving the pattern (Synchronous Update)

% define constants
N = 50; % size of network
T = 100; % simulation length

% generate random pattern (V) and weights matrix (W)
[V, W] = patternWeight(N);

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
%overlapCELL = cell(1, length(proportionNoise));

for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
    S = V; % re-initialize S (don't need this prob)
    x = addNoise(S, N, proportionNoise(noise_iter));
    h = zeros(1,N);
    storeH = cell(1,100);
    overlap =[];
    m = (x * V')/N; % initial overlap
    overlap = [overlap, m]

    for time = 2:T % iterate over each time point (why am I starting at 2?)
        for neuron = 1:N % iterate over each neuron in the network 
            % compute the *input potential* of a neuron for current neuron
            h(neuron) = (W(neuron,:) * x')/N;
%             x(neuron) = sign(h(neuron));
            % x_prime(neuron) = sign(h(neuron));
        end
        storeH{1,time} = h;
        % measure similarity between current state, x(t) and pattern, V.
        % x_prime = x;
        x=sign(h);
        figure()
        plot(x)
        %x = h;
        m = (x*V')/N;
        overlap = [overlap, m];
    end
    
    OL(:,noise_iter) = overlap;

    % plot m (should remain at 1)
    % since the initial network is all 1s, the dynamics of k
    % converged to a fixed point corresponding to the pattern (V) which
    % is most similar to the initial state. (?)
   
%     subplot(5,2,noise_iter)
%     plot(OL(:,noise_iter))
%     title(sprintf('%dth percent noise',proportionNoise(noise_iter)))
end