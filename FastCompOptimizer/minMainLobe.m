function minMainLobe = minMainLobe(pairs, intervals)

    % Returns minimum value of mainlobes
    % in the sum of the autocorrelation of the complementary pairs in "pairs"
    
    % The first and second row form a complementary pair, the second and 
    % third row form a complementary pair, and so on. We want the sum
    % of the autocorrelations for each pair to be like a delta function.
    % The central peak of this delta function is the mainlobe value.
    %pairs = reshape(pairs, [5*2, 8]);
    
    numPairs = size(pairs,1)/2;
    minMainLobe = -1;
    
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
    
    % The input is a matrix where each two rows is a complementary pair
    for i = 1:numPairs
        % Find max value of current mainlobe
        currPair = pairs((2*i-1):2*i,:); % Get next pair
        %currACF = xcorr(currPair(1,:)) + xcorr(currPair(2,:));
        %currMain = max(abs(currACF));
        currMain = currPair(1,:)*currPair(1,:)' + currPair(2,:)*currPair(2,:)';
        
        % Store smallest mainlobe value so far
        if (currMain < minMainLobe || minMainLobe == -1)
            minMainLobe = currMain;
        end
    end          

end

