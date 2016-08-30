function width = fullWidthHM(dataIn)
% Calculate full width half maximum
% This finds the width of the maximal peak in the input data
% (how many sequence elements are above 0.5 of the maximum value in the
% central peak)

% Locate the central peak
[maxVal,maxValInd] = max(dataIn);

% Find number of values above maximum to right of maximum
% in central peak
numAboveRight = 0;
while(dataIn(maxValInd+numAboveRight+1) > maxVal/2)
   numAboveRight = numAboveRight + 1; 
end

% Find number of values above maximum to left of maximum
% in central peak
numAboveLeft = 0;
while(dataIn(maxValInd-(numAboveLeft+1)) > maxVal/2)
   numAboveLeft = numAboveLeft + 1; 
end

width = numAboveRight + numAboveLeft + 1;