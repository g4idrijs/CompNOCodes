function bound = welchBound( codeLen, numCodesPar, transEvents )
% Returns the optimal Welch bound given the specified code length, the 
% number of codes transmitted in parallel, and the the number of transmit events

% The sum of the magnitudes squared of each code pair needs to be 1.

% This bound is a lower bound on the maximum of the non-mainlobe
% autocorrelation, and the maximum value in the sum of cross correlations
% between codes. (ex. xcorr(A1,B1) + xcorr(A2,B2))

N = codeLen;
M = numCodesPar;
K = transEvents;
bound = sqrt((M/K -1) / (M*(2*N-1)-1));

end

