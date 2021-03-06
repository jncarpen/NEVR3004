function [patternVec, weightMat] = patternWeight(N, selfConn)
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