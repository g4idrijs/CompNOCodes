function maxXcorr = maxXcorr( pairs, intervals )
    % Returns maximum sum of pairwise cross correlation
    
    % The input is a matrix where each two rows is a complementary pair
    
    % If first pair is A1, A2 and the second is B1,B2 then we are
    % interested in the maximum value of xcorr(A1,B1) + xcorr(A2,B2).
    % We find this for all pairs of pairs, and then choose the maximum.
    %pairs = reshape(pairs, [5*2, 8]);
    
    numPairs = size(pairs,1)/2;
    maxXcorr = -1;
    pairChoices = nchoosek(1:numPairs,2);  
    
    % If using NRI, add intervals of zeros between non-zero bits
    if length(intervals)
        x = pairs;
        x_int = zeros(size(x, 1), size(x, 2)+sum(intervals(1:size(x, 2)-1)));
        k = 2;
        x_int(:, 1) = x(:, 1);
        for i = 2:size(x, 2)
            insert = [zeros(size(x, 1), intervals(i-1)) x(:, i)];
            x_int(:, k:k+size(insert, 2)-1) = insert;
            k = k+size(insert, 2);
        end
        pairs = x_int;
    end
    
    for i = 1:size(pairChoices,1)
        % We find the sum of pairwise cross correlations
        % Ex. for pair A and pair B, we sum xcorr(A1,B1) + xcorr(A2,B2)
        currChoices = pairChoices(i,:);
        firstPair = pairs((2*currChoices(1)-1):2*currChoices(1),:);
        secPair = pairs((2*currChoices(2)-1):2*currChoices(2),:);
        
        currXcorr = xcorr(firstPair(1,:),secPair(1,:)) + xcorr(firstPair(2,:),secPair(2,:));
        currMax = max(abs(currXcorr));        

        % Store the largest cross cross correlation value
        if (currMax > maxXcorr || maxXcorr == -1)
            maxXcorr = currMax;
        end
    end


end

