% FLASH_COMPLEMENTARY_SIM 
% Simulate a flash image with complemenatry codes at many focal zones
% VERSION 1.0, August 10, 2016
%
% This version uses BFT to beamform images

%path('/home/tjh/git/zemp_lab/Matlab/bft',path);
%path('/home/tjh/git/fieldII',path);

%% Simulation properties

% Define codes used for excitation
codes = cell(0);
codes{1}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
codes{1}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
codes{1}.focus = [0 0 55/1000];

no_lines = 257;

% Transducer properties
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]

rx_fnum         = 2.1;   % Receive f-number - Set to 0 to turn off
                       % constant F-number reconstruction

%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
pht_pos = [0 0 20;
           0 0 30;
           0 0 40;
           0 0 50;
           0 0 60;
           0 0 70;
           0 0 80;] / 1000;         %  The position of the phantom
pht_amp = 20*ones(7,1);      %  The amplitude of the back-scatter

% Calculate minimum and maximum samples and times for phantom
Rmax = max(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) + 5/1000;
Rmin = min(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) - 5/1000;
if (Rmin < 0) Rmin = 0; end;
Tmin = 2*Rmin / c; Tmax = 2*Rmax / c;
Smin = floor(Tmin * fs); Smax = ceil(Tmax * fs);
max_code_length = max(cellfun(@(x) max(length(x.code), length(x.ccode)), codes));
Smin_c = Smin; Smax_c = Smax + max_code_length + 1000;

no_rf_samples = Smax - Smin + 1;
no_rf_samples_c = Smax_c - Smin_c + 1;

%  Initialize Field II and BFT
field_init(-1);
bft_init;

%  Set Field and BFT parameters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);

bft_no_lines(no_lines);

% Create transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

xdc = bft_linear_array(no_elements, width, kerf);

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

% Set the apodization
xdc_apodization(xmt, 0, ones(1,no_elements))
xdc_apodization(rcv, 0, ones(1,no_elements))


%% Simulate transmit and receive

% Define receive focusing - dynamic focusing on all lines
xdc_center_focus(xmt,[0 0 0])
xdc_center_focus(rcv,[0 0 0])
xdc_focus_times(rcv, 0, zeros(1,no_elements));

d_x_line = width*no_elements / (no_lines-1);
x_line = -(no_lines-1) / 2 * d_x_line;
z_points = linspace(Rmin, Rmax, 100);
no_points = length(z_points);
T = z_points/c*2;
for i = 1 : no_lines
    % Focus everywhere in the line
    bft_center_focus([x_line 0 0], i);
    bft_dynamic_focus(xdc, 0, 0, i); 
    
    % Constant F-Num reconstruction
    if (rx_fnum == 0)
        % Use all elements if fnum = 0
        bft_apodization(xdc, 0, ones(1, no_elements), i);
    else
        rx_ap           = z_points/rx_fnum;
        no_active_rx    = min(round(rx_ap/(width+kerf)), no_elements);

        mid_element = round((x_line+d_x_line*(no_lines-1)/2)/(width)+0.5);
        start_element = round(max(mid_element - no_active_rx/2 + 1, 1));
        end_element = round(min(mid_element+no_active_rx/2, no_elements));

        apo_vector_rx = zeros(no_points, no_elements);
        for j = 1:no_points
            apo_vector_rx(j, start_element(j):end_element(j)) = 1;
        end
    
        bft_apodization(xdc, T, apo_vector_rx, i);
    end
  
    x_line = x_line + d_x_line;
end


%% Simulate imaging with each code
% rf_data generated from imaging with each code and its associated focus
% zone is added together to simulate imaging with all codes. This is then 
% repeated for complementary versions.

% Rf data matrix holds non-complementary and complementary simulated
% RF data.
rf_data_m = zeros(no_rf_samples_c, no_elements, 2);
for i = 1:2*length(codes)
    if (i > length(codes))
        code = codes{i-length(codes)}.ccode;
        focus = codes{i-length(codes)}.focus;
        col = 2;
    else
        code = codes{i}.code;
        focus = codes{i}.focus;
        col = 1;
    end
    xdc_focus(xmt, 0, focus);
    xdc_excitation(xmt, code);
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);

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

    rf_data_m(1:size(rf_data, 1), :, col) = rf_data_m(1:size(rf_data, 1), :, col) + rf_data;
end

%% Decode images
% For each code, cross-correlate rf_data_sum with each code and add the
% result to rf_data_decoded
rf_data_nc_decoded = 0;
rf_data_c_decoded = 0;
for i = 1:length(codes)
    %temp_decoded = xcorr2(rf_data_m(:, :, 1), codes{i}.code');
    %temp_decoded = temp_decoded(length(codes{i}.code):end, :);
    temp_decoded = conv2(rf_data_m(:, :, 1), rot90(conj(codes{i}.code'), 2), 'valid');
    temp_decoded = temp_decoded(1:no_rf_samples, :);
    rf_data_nc_decoded = rf_data_nc_decoded + temp_decoded;

    %temp_decoded = xcorr2(rf_data_m(:, :, 2), codes{i}.ccode');
    %temp_decoded = temp_decoded(length(codes{i}.ccode):end, :);
    temp_decoded = conv2(rf_data_m(:, :, 2), rot90(conj(codes{i}.ccode'), 2), 'valid');
    temp_decoded = temp_decoded(1:no_rf_samples, :);
    rf_data_c_decoded = rf_data_c_decoded + temp_decoded;
end

% Can either add decoded data for non-complementary and complementary then
% beamform.
rf_data_decoded = rf_data_nc_decoded + rf_data_c_decoded;
bf_temp = bft_beamform(Tmin, rf_data_decoded);

% Or we can beamform decoded data then add results
%bf_temp = bft_beamform(Tmin, rf_data_nc_decoded);
%bf_temp = bf_temp + bft_beamform(Tmin, rf_data_c_decoded);

% Perform envelope detection and normalize
env_bf = abs(hilbert(bf_temp));
env_bf = env_bf / max(max(env_bf));

% Plot image
figure;
imagesc([-1/2 1/2]*no_lines*d_x_line*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);

