% Makes sure the data we get is from the right time
% to image the desired range between R_min and R_max

function rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements)

    start_sample = floor(start_time * fs)+1;
    end_sample = start_sample + size(rf_data, 1) - 1;

    if (start_sample < Smin_c)
      rf_data = rf_data(Smin_c-start_sample+1:end, :);
      start_sample = Smin_c;
    end
    rf_data = [zeros(start_sample - Smin_c, no_elements); rf_data; zeros(Smax_c - end_sample, no_elements)];

    if size(rf_data, 1) > no_rf_samples_c
      rf_data = rf_data(1:no_rf_samples_c, :);
    end