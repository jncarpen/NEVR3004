%% Check for latching dynamics
% Last modified: June 8, 2020

N = 50; % network size
P = 2; % number of patterns
alpha = P/N; % compute alpha value
timeSteps = 1500; % simulation length 
total_iterations = 1000; % number of complete simulations to run
allPatterns = zeros(P,N); % intialize cell array to store multiple patterns
W = zeros(N,N); % initalize empty weights matrix

for patt = 1:P
    [patternVec, weightMat] = patternWeight(N,0); % use function 'patternWeight' 
    allPatterns(patt,:) = patternVec; % store current pattern in matrix
    W = W + weightMat; % add weight matrices together 
end

V = allPatterns(1,:); % pull out pattern V
U = allPatterns(2,:); % pull out pattern U
pearsons = corrcoef(U,V); % compute correlation of the two patterns
convergence = cell(length(total_iterations));

for simulationIter = 1:total_iterations
    x = [V(1:25), U(26:50)]; % make cue half patV and half patU
    m(1) = (x * V')/N; % calculate overlap between state & pattern
    overlap = []; % initialize an overlap matrix
    m2overlap = [];
    
    for simLength = 1:timeSteps % total # of iterations *

        % choose between non-sequential or sequential update
        neuron = randi([1 50]); % random number from 1-50
        % neuron = repSeq(simLength); % use for asynchronous sequential update

        h(neuron) = (W(neuron,:) * x')/N; % compute input potential of neuron_i
        x(neuron) = sign(h(neuron)); % update the state of the network
        m(neuron) = (x*V')/N; % similarity between current state & pattern
        m2(neuron) = (x*U')/N; % similarity between current state & another pattern
        overlap = [overlap, m(neuron)]; %store in overlap matrix
        m2overlap = [m2overlap, m2(neuron)];
    end
    convergence{simulationIter} = [overlap]';
    m2convergence{simulationIter} = [m2overlap]';
end



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