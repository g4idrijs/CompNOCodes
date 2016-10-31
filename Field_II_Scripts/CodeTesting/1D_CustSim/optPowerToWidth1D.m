% Plot (peak power level)/(width above threshold) - call this clarity ratio
% Maximize how much the signal stands out above noise, and how sharp it is

% Threshold for determining central width
dBThresh = -45; % dB

% Fineness of simulation (numelements in impulse response)
len = 1000;

% Hold clarity ratio values
clarRatVect = [];

progCount = 1; % Keep track of progress on creating plot data

maxValNoise = 0; % Maximum value of uniform noise added

% Codes used (complementary pair)
codes = load('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Complementary Pairs\PR_CompPairs_VarLen\pairLen100.mat');
codes = codes.('x');
codes = codes./max(max(abs(codes)))*3; % Normalize the codes to [-3  3]

elRepeatRange = round(linspace(1, len, 100));
for custTimes = elRepeatRange;
    plotResult = 0; % Don't plot inside loop
    [highWidth, maxStrength] =  oneDimCustMain(custTimes, len, plotResult,dBThresh, maxValNoise, codes);
    clarRatVect = [clarRatVect maxStrength/highWidth];    
    
    if(mod(progCount,100) == 0) 
        disp(progCount);
    end
    progCount = progCount + 1;
end

%% Plot results
plot(elRepeatRange/len,clarRatVect,'-');
xlabel('Time per Element (time impulse response)')
ylabel('Peak power level / width above threshold')
title(strcat('Element Duration Effect on Clarity of Self-Decoding', ' Noise: ', num2str(maxValNoise)))

%% Get index of maximum
[maxVal, maxInd] = max(clarRatVect);
optTimePerEl = elRepeatRange(maxInd)/len; % In units of impulse response length
disp('Time per element for most clarity: (fraction of impulse response time)')
disp(optTimePerEl);