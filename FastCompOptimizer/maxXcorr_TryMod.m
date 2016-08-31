function maxXcorr = maxXcorr_TryMod( x, xc )
    % Returns maximum sum of pairwise cross correlation
    
    % The input is a matrix where each two rows is a complementary pair
    
    % If first pair is A1, A2 and the second is B1,B2 then we are
    % interested in the maximum value of xcorr(A1,B1) + xcorr(A2,B2).
    % We find this for all pairs of pairs, and then choose the maximum.
        
    pairs = xc;
    numPairs = size(pairs,1)/2;
    totCorr = zeros(1,2*size(pairs,2)-1);
    pairChoices = nchoosek(1:numPairs,2);  
    
    maxXcorr = 0;
    for i = 1:size(pairChoices,1)
        % We find the sum of pairwise cross correlations
        % Ex. for pair A and pair B, we sum xcorr(A1,B1) + xcorr(A2,B2)
        currChoices = pairChoices(i,:);

        firstPairRow = (2*currChoices(1)-1):2*currChoices(1);
        secPairRow = (2*currChoices(2)-1):2*currChoices(2);

        firstPair = pairs((2*currChoices(1)-1):2*currChoices(1),:);
        secPair = pairs((2*currChoices(2)-1):2*currChoices(2),:);
                
        currXcorr = xcorr(firstPair(1,:),secPair(1,:))+ xcorr(firstPair(2,:),secPair(2,:));             

        %totCorr = totCorr + currXcorr;
        maxXcorr = max(abs(currXcorr))+maxXcorr;
    end

    %maxXcorr = max(abs(totCorr));
end


