function [ metric, message ] = maxToAvg(sigData, clutterData)
    % Method of quantifying signal to clutter ratio
    % Returns max(signal)/avg(clutter)

    % Find max value in clutter region (dB)
    maxCent = max(max(sigData));

    % Find average value in clutter region (dB)
    avgClutter = mean(clutterData);
    
    % Take ratio of values
    metric = maxCent/avgClutter;

    message = 'Max Signal to Average Clutter Ratio:';
    
end


