% Makes sure the data we get is from the right time
% to image the desired range between R_min and R_max

function rf_data = alignRF_DavidMod(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements,printStuff,manualTrim)
        
    % start_time = the start time for the first sample (row) of rf_data
    
    % Find sample indices at start and end
    Ts = 1/fs; % Sampling period
    startSampInd = start_time/Ts+1; 
    endSampInd = startSampInd + size(rf_data, 1) - 1;
    
     rf_data = rf_data(manualTrim:end,:);     
    if(printStuff == 1) % On first transmit of a code
        disp('start_sample')
        disp(startSampInd)
        disp('end_sample')
        disp(endSampInd)
        
%         figure
%         imagesc(rf_data)
%         
%         title('Comparing Data Arrival Time')
    end       
    
    % We only want to collect data during a certain section of samples
    % Trim off extra data at the start   
%     if (startSampInd < Smin_c)
%       rf_data = rf_data(Smin_c-startSampInd+1:end, :);
%       startSampInd = Smin_c;
%     end
    
    % Add extra zeros at the start if our data doesn't start until later
    % than the standard.
    % Add extra zeroes at the end if our data ends before the standard
    % time.
    % rf_data = [zeros(startSampInd - Smin_c, no_elements); rf_data; zeros(Smax_c - endSampInd, no_elements)];
    
    % Just add zeroes at the end
    rf_data = [rf_data; zeros(no_rf_samples_c - size(rf_data,2), no_elements)];

    % Trim off extra data at the end if needed (never runs?)
    if size(rf_data, 1) > no_rf_samples_c
      rf_data = rf_data(1:no_rf_samples_c, :);
    end    
