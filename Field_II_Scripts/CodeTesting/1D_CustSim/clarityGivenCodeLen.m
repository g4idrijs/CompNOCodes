% Plot (peak power level)/(width above threshold) - call this clarity ratio
% Maximize how much the signal stands out above noise, and how sharp it is
% We transmit each element at the optimal amount of time: 0.25 of the
% impulse response length and measure clarity.

% Threshold for determining central width
dBThresh = -45; % dB

% Fineness of simulation (numelements in impulse response)
len = 1000;

% Hold clarity ratio values
clarRatVect = [];

progCount = 1; % Keep track of progress on creating plot data

maxValNoise = 2000; % Maximum value of uniform noise added

codeLenRange = 2:100;
for codeLen = codeLenRange
    
    % Number of times to try a code of length codeLen (and then average)
    numTimesPerLen = 50;
    clarRatSubVect = [];
    for codeRepInd = 1:numTimesPerLen
        codes = genCompPair_1D(ones(1,codeLen-1));
        codes = codes./max(max(abs(codes)))*3; % Normalize the codes to [-3  3]

        plotResult = 0; % Don't plot inside loop

        % Calculate clarity given current code
        custTimes = round(len/4); % Set optimal element repetition (for 1D)
        [highWidth, maxStrength] =  oneDimCustMain(custTimes, len, plotResult,dBThresh, maxValNoise, codes);    
        clarRatSubVect = [clarRatSubVect maxStrength/highWidth];
    end
    % Take average of results (sending different codes of same length)
    clarRatVect = [clarRatVect mean(clarRatSubVect)];      
    
    % Display progress
    disp(progCount);
    progCount = progCount + 1;
end
%% Plot results
plot(codeLenRange,clarRatVect,'-');
xlabel('Code length (number of elements in pattern)')
ylabel('Peak power level / width above threshold')
title(strcat('Code Length Effect on Clarity of Self-Decoding', ' Noise: ', num2str(maxValNoise)))

%% Get index of maximum
[maxVal, maxInd] = max(clarRatVect);
optTimePerEl = codeLenRange(maxInd); % In units of impulse response length
disp('Code length for most clarity: (fraction of impulse response time)')
disp(optTimePerEl);