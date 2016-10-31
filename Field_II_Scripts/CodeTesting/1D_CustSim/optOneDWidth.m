% Plot the width above dBThresh dB as a function of how long
% we transmit each element relative to impulse response length

% Threshold for determining central width
dBThresh = -45; % dB

% Fineness of simulation (numelements in impulse response)
len = 1000;

% Hold width above dBThresh dB as function of element transmit time
widthFunElemTime= [];

progCount = 1; % Keep track of progress on creating plot data

maxValNoise = 2000;

elRepeatRange = round(linspace(len/100, len*0.5, 1000));
for custTimes = elRepeatRange;
    plotResult = 0; % Don't plot inside loop
    widthFunElemTime = [widthFunElemTime oneDimCustMain(custTimes, len, plotResult,dBThresh,maxValNoise)];    
    
    disp(progCount);
    progCount = progCount + 1;
end

%% Plot results
plot(elRepeatRange/len,widthFunElemTime,'-');
xlabel('Time per Code Element (time impulse response)')
ylabel('Width Above -45 dB (time impulse response)')
title(strcat('Element Duration Effect on Sharpness of Self-Decoding', ' Noise: ', num2str(maxValNoise)))