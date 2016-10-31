function plotOneDim(F,titleStr, len)

% Normalize
norm = @(F)(F./max(max(abs(F))));

% Convert to dB
dB = @(F)20*log10(abs(F) + eps);

% Label x axis according to fraction of length of trasducer impulse
% response
xVect = (1:numel(F))/len;

% Plot
plotdB = 1;
if(plotdB == 1)
    toDisp = dB(norm(F));    
    yText ='Magnitude Output (dB)';
    plot(xVect, toDisp,'-')
    ylim([-50 0])
else
    toDisp = id(F);
    yText = 'Magnitude Output';
    plot(xVect, toDisp,'-')
end

ylabel(yText)

xlabel('Time (ImpResp length)')
title(titleStr)

