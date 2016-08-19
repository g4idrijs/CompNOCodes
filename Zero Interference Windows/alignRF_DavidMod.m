% Makes sure the data we get is from the right time
% to image the desired range between R_min and R_max

function rf_data = alignRF_DavidMod(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements)
    
    % start_time = the start time for the first row of rf_data
    
    % Figure out which samples we have in rf_data
    % (where the first sample was collected at time t = 1/fs)
    start_sample = floor(start_time * fs)+1;
    end_sample = start_sample + size(rf_data, 1) - 1;    
    
    % We only want to collect data during a certain section of samples
    % Trim off extra data at the start
    if (start_sample < Smin_c)
      rf_data = rf_data(Smin_c-start_sample+1:end, :);
      start_sample = Smin_c;
    end
    
    % Add extra zeros at the start if our data doesn't start until later
    % than the standard.
    % Add extra zeroes at the end if our data ends before the standard
    % time.
    rf_data = [zeros(start_sample - Smin_c, no_elements); rf_data; zeros(Smax_c - end_sample, no_elements)];

    % Trim off extra data at the end if needed (never runs?)
    if size(rf_data, 1) > no_rf_samples_c
      rf_data = rf_data(1:no_rf_samples_c, :);
    end    
