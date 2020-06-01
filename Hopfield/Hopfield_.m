%%%%%%%%%%%%%%%%%%%%%%%%
% HOPFIELD
%%%%%%%%%%%%%%%%%%%%%%%%

% Each neuron is a perceptron with +1/-1 output
% Each neuron *recieves input* from every other neuron 
% Each neuron *sends output* to every other neuron

% The network is symmetric where w_ij = w_ji
% A neuron "flips" if the weighted sum of other neuron's inputs is
% *opposite sign*

% Construct a default network and generate a pattern matrix
% initialize a network of N neurons (network size)
N = 50;

% initialize number of patterns
% nPatterns = 20;

% rnadomize an array of size N representing the
% binary state of each neuron (+1 or -1).
S = ones(1, 50);

% initialize a pattern, V (method #1)
V = rand(1, 50);
V = 2 * (V > 0.5) - 1; 

% initialize a random pattern of 0s and 1s (method #2)
% V = randi([0 1], 1,50); % initialize row vector of 1s and 0s
% V(V==0) = -1; % change all 0s to -1s

% initialize /synaptic/ weights matrix, W.
W = V' * V;

% set diagonal elements (self-connections) of W to 0.*
W = W - diag(diag(W));


% Task I: Establishing pattern stability

define the length of the simulation (# of updates).
T = 50;

set an initial time
t = 1;

initialize network by setting the state of each neuron
to the value it takes in the pattern V.
S = V;

calculate the overlap/similarity between state of the network, S
and the stored pattern(s), V.
m(1) = (S * V')/N;


compute field to each neuron that is h_i = sum_j(w_ij * s_j)
this is the weighted sum of outputs from all other neurons in network

for ti = 2:T % iterate over each time point (why am I starting at 2?)
    for neuron = 1:N % iterate over each neuron in the network 
    h(i) = sum((W(i,:) * S')/50); % 

        compute the *input potential* of a neuron for current time
        h(neuron) = (W(neuron,:) * S')/N;
        
        update the state of that neuron
        is opposite to the neuron's current state, S(neuron) will *flip*
        S(neuron) = sign(h(neuron));
        
    end
    
    measure similarity between current state, S(t) and pattern, V.
    m(ti) = (S * V')/N;
end

plot m (should remain at 1)
since the initial network is all 1s, the dynamics of the network
converged to a fixed point corresponding to the pattern (V) which
is most similar to the initial state. (?)

figure()
plot(m)

%  synchronous update i 
% define constants
N = 50; % size of network
T = 50; % simulation length
t = 1; % initial time

% generate random pattern (V) and weights matrix (W)
[V, W] = patternWeight(N);

% initialize network by setting the state of each neuron
% to the value it takes in the pattern V.
%S = V;
% 
% % randomly add noise to the state matrix
% x = addNoise(S, N, 0.6); % function below
% 
% % calculate the overlap/similarity between state of the network, S
% % and the stored pattern(s), V.
% m(1) = (x * V')/N;

% compute field to each neuron that is h_i = sum_j(w_ij * s_j)
% this is the weighted sum of outputs from all other neurons in network

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];

for noise_iter = proportionNoise % loop through all possible proportions of noise
    S = V; % re-initialize S (don't need this prob)
    x = addNoise(S, N, noise_iter);
    overlap = [];
    h = zeros(1,N);
    storeH = cell(1,100);
    % calculate the overlap/similarity between state of the network, S
    % and the stored pattern(s), V.
    m(1) = (x * V')/N;
    % non-sequential (random) synchronous update
    for time = 2:100 % iterate over each time point (why am I starting at 2?)
        for neuron = 1:N % iterate over each neuron in the network 
            % compute the *input potential* of a neuron for current neuron
            h(neuron) = (W(neuron,:) * x')/N;
%             x(neuron) = sign(h(neuron));
            % x_prime(neuron) = sign(h(neuron));
        end
        storeH{1,time} = h;
        % measure similarity between current state, x(t) and pattern, V.
        % x_prime = x;
        % x=sign(h);
        m(time) = (x*V')/N;
        overlap = [overlap, m(time)];
    end

    % plot m (should remain at 1)
    % since the initial network is all 1s, the dynamics of the network
    % converged to a fixed point corresponding to the pattern (V) which
    % is most similar to the initial state. (?)
   
    noiseIdx = find(proportionNoise == noise_iter);
    subplot(5,2,noiseIdx)
    plot(m)
    title(sprintf('%dth percent noise',noise_iter*100))
end

% asynchronous non-sequential update

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


% asynchronous sequential update

N = 50; % size of network
[V, W] = patternWeight(N); % generate random pattern & weights matrix
iterations = 200; % simulation length
repSeq = repmat(1:N, 1,(iterations/N)); % creating a repeating matrix

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
figure() % initialize figure

for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
    S = V; % initialize network
    x = addNoise(S, N, proportionNoise(noise_iter)); % create vector w/ varied amounts of noise
    m(1) = (x * V')/N; % calculate overlap between state & pattern
    overlap = []; % initialize an overlap matrix
       
        for simNum = 1:iterations % total # of iterations *
            neuron = repSeq(simNum);
            h(neuron) = (W(neuron,:) * x')/N; % compute input potential of neuron_i
            x(neuron) = sign(h(neuron)); % update the state of the network
            m(neuron) = (x*V')/N; % similarity between current state & pattern
            overlap = [overlap, m(neuron)]; %store in overlap matrix
        end
        
        subplot(5,2,noise_iter)
        plot(overlap)
        title(sprintf('%dth percent noise',proportionNoise(noise_iter)))
end

%% synchronous update

N = 50; % size of network
[V, W] = patternWeight(N); % generate random pattern & weights matrix
iterations = 200; % simulation length

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
figure() % initialize figure

for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
    S = V; % initialize network
    overlap = []; % initialize an overlap matrix
    H = zeros(iterations,N); % input potential
    X = zeros(iterations,N); % state of network
    %X(1,:) = addNoise(S, N, proportionNoise(noise_iter)); % create vector w/ varied amounts of noise
    %m(1) = (X(1,:) * V')/N; % calculate overlap between state & pattern

    
        for simNum = 1:iterations % total # of iterations *
            if simNum == 1
                X(simNum,:) = addNoise(S, N, proportionNoise(noise_iter)); % create vector w/ varied amounts of noise
                m(simNum) = (X(simNum,:) * V')/N; % calculate overlap between state & pattern
            else
            for neuron = 1:N
            H(simNum, neuron) = (W(neuron,:) * X(simNum-1,:)')/N; % compute input potential of neuron_i
            end
            X(simNum,:) = sign(H(simNum,:)); % update the state of the network
            m(simNum) = (X(simNum,:)*V')/N; % similarity between current state & pattern
            overlap = [overlap, m(simNum)]; %store in overlap matrix
            end
        end
        
        subplot(5,2,noise_iter)
        plot(overlap)
        title(sprintf('%dth percent noise',proportionNoise(noise_iter)))
end

% synchronous update ii

N = 50; % size of network
[V, W] = patternWeight(N); % generate random pattern & weights matrix
iterations = 200; % simulation length

% create a vector of proportions
proportionNoise = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
figure() % initialize figure

for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
    S = V; % initialize network
    overlap = []; % initialize an overlap matrix
    H = cell(1,iterations); % input potential
    X = cell(1,iterations); % state of network
    m = zeros(1,iterations); % initialize an overlap matrix
    
        for simNum = 1:iterations % total # of iterations *
            if simNum == 1
                X{simNum} = addNoise(S, N, proportionNoise(noise_iter))'; % create vector w/ varied amounts of noise
                m(simNum) = (V * X{simNum})./N; % calculate overlap between state & pattern
                overlap = [overlap m(simNum)];
            else
            h = (W*X{simNum-1}); % compute input potential of neuron_i
            H{simNum} = h./N;
            X{simNum} = sign(H{simNum}); % update the state of the network
            m(1,simNum) = (V * X{simNum})./N; % similarity between current state & pattern
            overlap = [overlap, m(simNum)]; %store in overlap matrix
            end
        end
        
        subplot(5,2,noise_iter)
        plot(overlap)
        title(sprintf('%dth percent noise',proportionNoise(noise_iter)))
end
% store multiple patterns

generate a bunch of patterns and corresponding weight matrices
N = 50;
number_patterns = 1;
allPatterns = zeros(number_patterns,N);
W = zeros(N,N);

for patt = 1:number_patterns 
    [patternVec, weightMat] = patternWeight(N,3);
    allPatterns(patt,:) = patternVec;
    W = W + weightMat; % add weight matrices together
end

% asynchronous non-sequential update (test w/ multiple patterns)

create a vector of proportions
proportionNoise = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1]; 
total_iterations = 1000;
figure() % initialize figure
V = allPatterns(1,:);
convergence = cell(length(total_iterations),length(proportionNoise));
    
for simulationIter = 1:total_iterations
    for noise_iter = 1:length(proportionNoise) % loop through all possible proportions of noise
        S = V; % initialize network
        x = addNoise(S, N, proportionNoise(noise_iter)); % create vector w/ varied amounts of noise
        m(1) = (x * V')/N; % calculate overlap between state & pattern
        overlap = []; % initialize an overlap matrix


            for simLength = 1:500 % total # of iterations *

                neuron = randi([1 50]); % random number from 1-50
                h(neuron) = (W(neuron,:) * x')/N; % compute input potential of neuron_i
                x(neuron) = sign(h(neuron)); % update the state of the network
                m(neuron) = (x*V')/N; % similarity between current state & pattern
                overlap = [overlap, m(neuron)]; %store in overlap matrix
            end

            convergence{simulationIter, noise_iter} = [overlap]';

            subplot(3,4,noise_iter)
            plot(overlap)
            title(sprintf('%d percent noise',proportionNoise(noise_iter)*100))
    end
end

% convergence

for i = 1:total_iterations
    for j = 1:length(proportionNoise)
        minIter = find(convergence{i,j}==1, 1); % find the min iteration that the network converged
        if isempty(minIter) % if it never converges
            minIter = NaN; % set convergence value to NaN
        end
        minCon(i,j) = minIter; % store it
    end
end

% plot convergence histograms

figure()
for col=1:11
    subplot(3,4,col)
    histogram(minCon(:,col), 10)
end

       
%
%%%%%%%%%%%%%%%%%%%
   FUNCTIONS
%%%%%%%%%%%%%%%%%%%

generate random noise

function noisyState = addNoise(stateVector, N, perNoise)
    
    FUNCTION DESCRIPTION
    This function will take in a 1xN stateVector, where N is the size of the
    network, a scalar value perNoise which is the proportion of noise to
    add to the stateVector. The function will output a 1xN noisyState
    vector. If perNoise=0, noisyState will just be the original
    stateVector.
    
    if perNoise == 0
        noisyState = stateVector;
    else    
    permute = randperm(N, round(N*perNoise));
    for j = permute
        R = rand();
        if R > 0.5
            R=1;
        else
            R=-1;
        end
    stateVector(j) = R;
    noisyState = stateVector;
    end
    end
end

% generate pattern and weight matrices

function [patternVec, weightMat] = patternWeight(N, selfConn)
generate pattern and weight matrices

initialize a pattern
patternVec = rand(1,N);
patternVec = 2*(patternVec>0.5)-1;

initialize /synaptic/ weights matrix, W (divide by N?)
weightMat = (patternVec' * patternVec);

if selfConn == 0
    set diagonal elements (self-connections) of W to 0.*
    weightMat = weightMat - diag(diag(weightMat));
elseif selfConn == 1
    set self-connections = 1
    weightMat(logical(eye(size(weightMat)))) = 1; 
elseif selfConn == 2
    set self-connections = -1
    weightMat(logical(eye(size(weightMat)))) = -1;
else
end
end