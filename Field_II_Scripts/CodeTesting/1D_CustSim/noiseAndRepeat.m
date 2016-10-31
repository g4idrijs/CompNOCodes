len = 6000;  % Simulation fineness
maxValNoise = 2000*len;
dBThres = -45; % dB (threshold for center width)

% repFrac of 0.255 maximizes clarity

% figure
% figCount = 1;
% for repFrac = linspace(0.1, 0.5, 8) % Fraction of impulse response to repeat code
%     subplot(2,4,figCount)
%     [highWidth, maxStrength]  = oneDimCustMain(round(repFrac*len), len, 1, dBThres, maxValNoise)    
%     
%     title(sprintf(' elemRepeat: %.2f. width: %.2f.', repFrac, highWidth));
%     
%     figCount = figCount + 1;
% end

% Illustrate what happens when we don't repeat
figure
figCount = 1;
for maxValNoise = exp(linspace(1,7.5,8))
    subplot(2,4,figCount)
    [highWidth, maxStrength]  = oneDimCustMain(1, len, 1, dBThres, maxValNoise);  
    
    title(sprintf('Noise: %.2f', maxValNoise));
    
    figCount = figCount + 1;
end