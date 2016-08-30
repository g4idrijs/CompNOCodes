function tagMat = tagMatrix( intVect )
% Generates a matrix with "1" where we have a tag bit and "0" elsewhere
% intVect tells us the interval between tag bits
% Ex. tagMatrix([1,2,3]) should generate:
% [1101001
% 1101001
% 1101001
% 1101001]

% We add n-1 zeroes for a value of n in the intVect vector.
% We have lenght(intVect) + 1 ones
lenTagVect = sum(intVect - 1) + length(intVect) + 1;

tagVect = zeros(1,lenTagVect);

% Always start with a 1
tagVect(1) = 1;

% Use cumulative sums to know where to put the next 1
tagVect(1+cumsum(intVect)) = 1;

% Create matrix of appropriate size
% Each row is a copy of tagVect
tagMat = repmat(tagVect, [length(intVect) + 1,1]);

end

