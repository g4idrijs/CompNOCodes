function [ metric, message ] = maxToMax( sigData, clutterData)
    % Method of quantifying signal to clutter ratio
    % Returns max(signal)/max(clutter)

    % Find max value in clutter region (dB)
    maxCent = max(max(sigData));

    % Find max value in clutter region (dB)
    maxClutter = max(max(clutterData));
    
    % Take ratio of values
    metric = maxCent/maxClutter;

    message = 'Max Signal to Max Clutter Ratio';
    
end

