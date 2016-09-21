clear

addpath('C:\Users\User\Dropbox\Grad_School\HighSpeedCodes\CompNOCodes\FastCompOptimizer')

%% Load in the appropriate pairs
numPairs = 16;
N = 10;

% tempDat = load('compPairs_len_10_simMain.mat');
% x = tempDat.('pairsSoFar');

tempDat = load('len10_16codes_minInt2');
x= tempDat.('x');



%% Compare to Welch bound
% disp('Welch bound:')    
% disp(welchBound( N, numPairs, 2 ))

% Would like to be welch bound
% Max cross correlation of normalized codes
% We set the magnitude of each individual pair to sqrt(2)
disp('Welch metric:') 
intervals = [];
disp(maxXcorr(normr(x)/sqrt(2),intervals))


%% Plot autocorrelation
% Create x axis (shift from -(codeLength-1) to +(codeLength-1)
xAxis = -(N-1):(N-1);
figure
minACFMain = -1;
for i=1:numPairs
    currPair = x((2*i-1):2*i,:);
    currAc = xcorr(currPair(1,:)) + xcorr(currPair(2,:));
    [maxVal, whereMax] = max(currAc);
    currAc(whereMax) = 0;
    sideMax = max(currAc);
    
    % mainToSide = maxVal/sideMax; % Check main-sidelobe in ACF is good
    
    plot(xAxis,xcorr(currPair(1,:)) + xcorr(currPair(2,:)));        
    hold on
    
    % Keep track of smallest mainlobe value
    if (maxVal < minACFMain || minACFMain == -1)
        minACFMain = maxVal;
    end
end    

% Min ACF mainlobe to CCF sidelobe ratio
ACFToCCF = minACFMain / maxXcorr(x, intervals);
disp('ACF to CCF:');
disp(ACFToCCF);

% Plot pairwise sum of cross correlation
plotCcSum = 1;

if (plotCcSum == 1 && numPairs > 1)

    pairChoices = nchoosek(1:numPairs,2);
    for i = 1:size(pairChoices,1)
        % We find the sum of pairwise cross correlations
        % Ex. for pair A and pair B, we sum xcorr(A1,B1) + xcorr(A2,B2)
        currChoices = pairChoices(i,:);
        firstPair = x((2*currChoices(1)-1):2*currChoices(1),:);
        secPair = x((2*currChoices(2)-1):2*currChoices(2),:);

        currXcorr = xcorr(firstPair(1,:),secPair(1,:)) + xcorr(firstPair(2,:),secPair(2,:));
        plot(xAxis,currXcorr)
        hold on
    end
    title(sprintf('Optimized: ACF and CCF with Code Length = %d and #Pairs = %d', N,numPairs));
else      
    title(sprintf('ACF with Code Length = %d and #Pairs = %d', N,numPairs));
end

xlabel('Shift amount')
ylabel('Magnitude')    

% Record the ACF to CCF ratio
str = sprintf('ACF to CCF ratio: %1.2f',ACFToCCF);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
