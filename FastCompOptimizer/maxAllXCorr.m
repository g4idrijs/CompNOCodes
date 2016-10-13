% Consider all pairs of codes in pairs - not just complementary pairs
% Gets max_pairs(max_n(abs(xcorr(pair)))

% The input is a matrix where each two rows is a complementary pair
function maxAllXCorrVal = maxAllXCorr(pairs)

numCodes = size(pairs,1);
pairChoices = nchoosek(1:numCodes,2);  

% Maximum cross correlation seen so far
maxXcorr = -1;

for currPair = 1:size(pairChoices,1)
    currChoices = pairChoices(currPair,:);
    firstCode = pairs(currChoices(1),:);
    secCode = pairs(currChoices(2),:);

    currXcorr = xcorr(firstCode, secCode);
    currMax = max(abs(currXcorr));        

    % Store the largest cross cross correlation value
    if (currMax > maxXcorr || maxXcorr == -1)
        maxXcorr = currMax;
    end
end

maxAllXCorrVal = maxXcorr;