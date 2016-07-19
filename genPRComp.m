% Goal is to generate kind of random codes that
% have good cross correlation properties

merit = 0;
meritThresh = 4;

numPairs = 4;
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
    len = size(pairs,2);
end
a = 3;

%% (NOT MUCH DONE HERE YET)
% Use straightforward approach to generate small ones
% Either generalize this approach or use papers to make longer codes

% We need the autocorrelations to sum to zero at each shift (except center)
% This puts N-1 constraints on two codes of length N each

% So we have 2N - (N -1) = N + 1 free parameters
% Ex. with length 4 codes, we have 5 free parameters

% What if we set the first code entirely? (h1, h2, h3, h4)
% And we also set the first value of the second code: h1'.

% First equation F(h1,h4,h1',h4') -> solve -> h4'

% The second equation puts a constrain on h3' and h4'
% Then h2', h3', h4'
% Sort of like backpropogation idea

% Set first code randomly
codeMean = 0;
codeStd = 0.5;
codeLen = 4;
code1 = random('Normal',codeMean,codeStd,1,4);

% Set second code to make it complementary
