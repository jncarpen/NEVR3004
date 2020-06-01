
function noisyState = addNoise(stateVector, N, perNoise)
    % perNoise = 0.2 for example (20% noise)
    permute = randperm(N, round(N*perNoise));
    for j = permute
        R = rand();
        if R > 0.5;
            R=1;
        else
            R=-1;
        end
    stateVector(j) = R;
    noisyState = stateVector;
    end
end