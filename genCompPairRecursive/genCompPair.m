function [ finalPair ] = genCompPair( shiftList )
% Trying out recursive algorithm in "New Multilevel Complementary Pairs of Sequences"
% Sounds promising because it has arbitary real paramters we can set
% randomly

% "No overlapping or gaps" if permuted 0 2^[1..N]
% S = [0 2.^[3 2 1]];     % Shifting values (seems to set zero patterns)
S = [0 shiftList];
maxCodeLen = 1+sum(S);  % The longest length of codes defined using S

% Real numbers - hopefully can use to determine frequency properties
% A = [-1 2 3 -2 -1 1];  
seedMean = 0;
seedStd = 1;
A = random('Normal',seedMean,seedStd,1,length(S));

% Initial setup
numIter = length(S);
aPrev = [1 zeros(1,maxCodeLen - 1)];
bPrev = [0 zeros(1,maxCodeLen - 1)];
aNext = aPrev; bNext = bPrev;

for n = 1:numIter  
    % Get a copy of b shifted to the right
    shiftedB = [zeros(1,S(n)) bPrev(1:(length(bPrev) - S(n)))];
    
    % Create a longer complementary pair
    aNext = aPrev + A(n)*shiftedB;
    bNext = A(n)*aPrev-shiftedB;
 
    aPrev = aNext;
    bPrev = bNext;
end

finalPair = [aNext; bNext];
end

