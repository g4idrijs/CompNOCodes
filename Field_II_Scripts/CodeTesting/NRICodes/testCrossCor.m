% Idea:
% -Keep codes basically same (Hammard matrix generated)
% -However, allow up to two tag bits to line up at the same time
% -For example, 2,3,5 would now be allowed

% Generates the code matrix
% (Each row is a matrix)
% Note that the number of elements in the interval spacing vector
% must be one less than a multiple of 4 (to get Hadamard matrix)
% (The number of tag bits per row is one more than the number of intervals
% and the number of tag bits per row is the side size of the Hadamard
% matrix)
size = 7;
seed = [2,3,4];
sideLobeBound = 3;
intSeq = createIntVect(size,seed,sideLobeBound-1,[1])
V = hadamardCodes(intSeq)

% Calculate maximum cross correlation (ignoring central value for
% autocorrelation)
% The autocorrelation values are along the main diagonal
% It's interesting to note that no cancelling occurs in at least once case
% -Because one code is all 1's, so no cancellation occurs
crossCorrMat(V)







