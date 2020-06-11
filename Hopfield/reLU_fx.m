function [output] = reLU_fx(X)
%RELU_FX Summary of this function goes here
%   Detailed explanation goes here
    for i=X(1):length(X)
    if X >1
        output(i) =1;
    elseif X<0
        output(i)=0;
    end
end

