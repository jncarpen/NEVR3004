%% Asynchronous update:
% sequential & non-sequential, multiple patterns, variable
% self-connections, multiple iterations...

% Jo Carpenter
% Last modified: June 6, 2020

%% store multiple patterns

P = 12; % number of patterns to loop through
N = 50; % number of nodes in the network (total size)
convergence2 = cell(1,P); % create convergence array for outer loop
timeSteps = 500; % number of iterations in the main simulation
repSeq = repmat(1:N, 1,(timeSteps/N)); % create repeating matrix

for NP = 1:2:P % loop through every other pattern
    % generate a bunch of patterns and corresponding weight matrices
    number_patterns = NP;
    alpha(NP) = number_patterns/N; % calculate alpha value
    allPatterns = zeros(number_patterns,N); % matrix of all patterns generated
    W = zeros(N,N); % initialize empty weights matrix

    for patt = 1:number_patterns 
        
        % to change nature of self-connections see *patternWeight* function
        % at end of script ***
        
        [patternVec, weightMat] = patternWeight(N,2); % use function 'patternWeight' 
        allPatterns(patt,:) = patternVec; % store current pattern in matrix
        W = W + weightMat; % add weight matrices together
        
    end

    %% asynchronous update (test w/ multiple patterns)

    % create a vector of proportions
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

                for simLength = 1:timeSteps % total # of iterations *

                    % choose between non-sequential or sequential update
                    % rules (just comment out for now)
                    
                    neuron = randi([1 50]); % random number from 1-50
                    % neuron = repSeq(simLength); % use for asynchronous sequential update
                 
                    h(neuron) = (W(neuron,:) * x')/N; % compute input potential of neuron_i
                    x(neuron) = sign(h(neuron)); % update the state of the network
                    m(neuron) = (x*V')/N; % similarity between current state & pattern
                    overlap = [overlap, m(neuron)]; %store in overlap matrix
                end

                convergence{simulationIter, noise_iter} = [overlap]';
                
%                 subplot(3,4,noise_iter) % generate subplot for each noise iteration
%                 plot(overlap)
%                 title(sprintf('%d percent noise',proportionNoise(noise_iter)*100))

        end
    end
    
    convergence2{1,NP} = convergence; % save convergence values for each number of patterns
    
end

%% convergence

minCon2 = cell(1,P);
propCon2 = cell(1,P);

for NP=1:2:P
    for i = 1:total_iterations
        for j = 1:length(proportionNoise)
            minIter = find(convergence2{1,NP}{i,j}==1, 1); % find the min iteration that the network converged
            if isempty(minIter) % if it never converges
                minIter = NaN; % set convergence value to NaN
            end
            minCon(i,j) = minIter; % store it
        end
    end
    
        % determine proportion of times that the network converged for different levels of noise
        for col = 1:length(minCon(1,:))
            propConverged(col) = 1 - (sum(isnan(minCon(:,col)))/length(minCon(:,col)));
        end
        
        % save values for each pattern iteration
        minCon2{1,NP} = minCon; 
        propCon2{1,NP} = propConverged;
        
end


%% plot convergence histograms

% figure()
% for col=1:11
%     subplot(3,4,col)
%     histogram(minCon(:,col))
% end

%% remove gaps from files
minConData = {minCon2{1,1}, minCon2{1,3}, minCon2{1,5}, minCon2{1,7}, minCon2{1,9}, minCon2{1,11}};
ConvergenceData = {convergence2{1,1}, convergence2{1,3}, convergence2{1,5}, convergence2{1,7}, convergence2{1,9}, convergence2{1,11}};



%% FUNCTIONS:

%% 1. noisyState: generate random noise

function noisyState = addNoise(stateVector, N, perNoise)
    
    % FUNCTION DESCRIPTION
%     This function will take in a 1xN stateVector, where N is the size of the
%     network, a scalar value perNoise which is the proportion of noise to
%     add to the stateVector. The function will output a 1xN noisyState
%     vector. If perNoise=0, noisyState will just be the original
%     stateVector.
    
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

%% 2. patternWeight: generate pattern and weight matrices

function [patternVec, weightMat] = patternWeight(N, selfConn)
% FUNCTION DESCRIPTION
% generate pattern and weight matrices

% initialize a pattern
patternVec = rand(1,N);
patternVec = 2*(patternVec>0.5)-1;

% initialize /synaptic/ weights matrix, W (divide by N?)
weightMat = (patternVec' * patternVec);

if selfConn == 0
%     set diagonal elements (self-connections) of W to 0.*
    weightMat = weightMat - diag(diag(weightMat));
elseif selfConn == 1
%     set self-connections = 1
    weightMat(logical(eye(size(weightMat)))) = 1; 
elseif selfConn == 2
%     set self-connections = -1
    weightMat(logical(eye(size(weightMat)))) = -1;
else
    print("Error: selfConn specification must be between 0 and 2")
end
end