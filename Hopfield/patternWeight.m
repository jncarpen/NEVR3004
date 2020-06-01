%% generate pattern and weight matrices

function [patternVec, weightMat] = patternWeight(N)
%generate pattern and weight matrices

% initialize a pattern
patternVec = rand(1,N);
patternVec = 2*(patternVec>0.5)-1;

% initialize /synaptic/ weights matrix, W (divide by N?)
weightMat = (patternVec' * patternVec);
% set diagonal elements (self-connections) of W to 0.*
weightMat = weightMat - diag(diag(weightMat));

end
