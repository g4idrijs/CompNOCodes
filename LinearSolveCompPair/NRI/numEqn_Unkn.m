% Determine the number of equations and unknowns
% when trying to solve linearly for complementary codes
% using NRI codes

% Construct NRI code with desired zero interval lengths
 zeroInts = [1,2,3]; % Length of zero intervals, in order
 v = constructNRI(zeroInts);
 
% The number of equations will be the number of non-zero overlaps
% not including the central overlap
% (at each such non-zero overlap, we have an equation)
% (ex. x_1*x_10 + xc_1*xc_10 = 0)
N = numel(v); % Length of our code
numEqn = round(xcorr(v));
numEqn = nnz(numEqn(1:(N-1)))
 
% Find number of unknowns
numUnk = nnz(v)*2 % We double the number of unknowns, because we have two codes

% I feel like linearizing these equations will be impossible by 
% my previous proof