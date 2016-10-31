clear
figure

plotLog = 1;

% Normalize log plot
normLog = 1;

addpath('C:\Users\Zemp-Lab\Desktop\OvernightPairGeneration\GitCodes\CompNOCodes\FastCompOptimizer')
addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\FastCompOptimizer')

%% Load in the appropriate pairs

% tempDat = load('bestPairsSQP41');
tempDat = load('lowAllCC_2pairs_length100_15p8');

x = tempDat.('x');
numPairs = size(x,1)/2;
N = size(x,2);

%% Plot autocorrelation
% Create x axis (shift from -(codeLength-1) to +(codeLength-1)
xAxis = -(N-1):(N-1);

maxACF = -1;
for i = 1:numPairs
    % Keep track of largest mainlobe value
    currPair = x((2*i-1):2*i,:);
    currAc = xcorr(currPair(1,:)) + xcorr(currPair(2,:));
    
    maxVal =  max(abs(currAc));
    
    if(maxVal > maxACF)
        maxACF = maxVal;
    end
end


minACFMain = -1;
for i=1:numPairs
    currPair = x((2*i-1):2*i,:);
    currAc = xcorr(currPair(1,:)) + xcorr(currPair(2,:));
    [maxVal, whereMax] = max(currAc);
    currAc(whereMax) = 0;
    sideMax = max(currAc);
    
    % mainToSide = maxVal/sideMax; % Check main-sidelobe in ACF is good
    
    if(plotLog == 1)
        if(normLog == 1)
             plot(xAxis,20*log10(abs(xcorr(currPair(1,:)) + xcorr(currPair(2,:)))/maxACF + eps));  
        else
             plot(xAxis,20*log10(abs(xcorr(currPair(1,:)) + xcorr(currPair(2,:))) + eps));  
        end
    else
       plot(xAxis,xcorr(currPair(1,:)) + xcorr(currPair(2,:)));  
    end
           
    hold on
    
    % Keep track of smallest mainlobe value
    if (maxVal < minACFMain || minACFMain == -1)
        minACFMain = maxVal;
    end   
end    

% Min ACF mainlobe to CCF sidelobe ratio
ACFToCCF = minACFMain / maxXcorr(x,[]);
% ACFToCCF = minACFMain / maxXcorr(x, intervals);
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
        if(plotLog == 1)
            if(normLog == 1)
                plot(xAxis,20*log10(abs(currXcorr)/maxACF+eps));
            else
                plot(xAxis,20*log10(abs(currXcorr)+eps));
            end
        else
            plot(xAxis,currXcorr)
        end
        hold on
    end
    title(sprintf('ACF and CCF with Code Length = %d and #Pairs = %d', N,numPairs));
else      
    title(sprintf('ACF with Code Length = %d and #Pairs = %d', N,numPairs));
end

xlabel('Shift amount')
if(plotLog == 1)
   ylabel('Magnitude (dB)') 
else
   ylabel('Magnitude')  
end

% Record the ACF to CCF ratio
str = sprintf('ACF to CCF ratio: %1.2f',ACFToCCF);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
