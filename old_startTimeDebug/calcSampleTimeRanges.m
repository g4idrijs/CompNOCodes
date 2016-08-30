function [Rmax, Rmin, Tmin, Smin, max_code_length, Smin_c, Smax_c, no_rf_samples, no_rf_samples_c] =...
    calcSampleTimeRanges(pht_pos, codes,c,fs)

% Minimum and maximum distance to image
Rmax = max(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) + 5/1000;
Rmin = min(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) - 5/1000;
if (Rmin < 0) Rmin = 0; end;

% Minimum and maximum time to wait while imaging
Tmin = 2*Rmin / c; Tmax = 2*Rmax / c;

% Minimum and maximum samples to take
% (if we started collecting at time t = 0)
Smin = floor(Tmin * fs); Smax = ceil(Tmax * fs);
max_code_length = max(cellfun(@(x) max(length(x.code), length(x.ccode)), codes));

% we have to collect a few more sampls to get the code as it comes back
Smin_c = Smin; Smax_c = Smax + max_code_length + 1000;

% Total number of samples to take
no_rf_samples = Smax - Smin + 1;
no_rf_samples_c = Smax_c - Smin_c + 1;