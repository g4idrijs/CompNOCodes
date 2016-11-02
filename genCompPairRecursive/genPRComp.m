% Goal is to generate kind of random codes that
% have good cross correlation properties

merit = 0;
meritThresh = 4;

numPairs = 2;
while(merit < meritThresh)

    % Generate a few codes and see how they cross correlate
    shiftList = 1:40;
    pairs = [];
    
    for i = 1:numPairs
        pairs = [pairs ; genCompPair(shiftList)];
    end

    % Check complementarity
    % complResult = areCompl(pairs(1:2,:))
    % complResult = areCompl(pairs(3:4,:))

    % Find minimum main lobe value
    % This is the minimum main lobe value in the sum
    % of the ACFs of each code set
    numPairs = size(pairs,1)/2;
    minMainLobe = -1;
    for i = 1:numPairs
        currPair = pairs((2*i-1):2*i,:);
        currACF = xcorr(currPair(1,:)) + xcorr(currPair(2,:));
        currMain = max(abs(currACF));

        if (currMain < minMainLobe || minMainLobe == -1)
            minMainLobe = currMain;
        end
    end


    % Find the cross correlation sum values
    % There will be one of these for each pair of codes
    maxXcorr = -1;
    pairChoices = nchoosek(1:numPairs,2);
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
   
    merit = minMainLobe/maxXcorr;
    if (merit > 3)
        disp(merit);
    end
    
end
len = size(pairs,2)
