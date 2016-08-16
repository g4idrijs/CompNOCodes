function [beam, maxDelay] = focusBeam(excitation, focus, xElements, fs, c)
    delays = (focus(3)-sqrt((xElements-focus(1)).^2+focus(3)^2))./c * fs;
    delays = delays-min(delays);
    maxDelay = max(delays)/fs;
    %delays = delays-min(delays);
    %delays = delays-max(delays);
    %delays = delays-min(delays);
    
    maxLength = ceil(max(abs(delays))+length(excitation));
    beam = zeros(maxLength, length(delays));
    for i = 1:length(delays)
        %nPreZeros = round(delays(i));
        %nPostZeros = maxLength-nPreZeros-length(excitation);
        %beam(:, i) = [zeros(1, nPreZeros) excitation zeros(1, nPostZeros)];
        
        %beam{i} = delayseq(excitation, delays(i), fs);
        tmpBeam = interp1(1:length(excitation), excitation, (1-delays(i)):length(excitation), 'linear', 0);
        beam(1:length(tmpBeam), i) = tmpBeam;
        %beam{i} = tmpBeam;
    end
end