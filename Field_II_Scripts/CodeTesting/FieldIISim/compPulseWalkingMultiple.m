% FLASH_COMPLEMENTARY_SIM 
% Simulate a flash image with complemenatry codes at many focal zones
% VERSION 1.0, August 10, 2016
%
% This version uses BFT to beamform images

%path('/home/tjh/git/zemp_lab/Matlab/bft',path);
%path('/home/tjh/git/fieldII',path);
tempLoad = load(['../../Complementary Pairs/compPairs_len_10_simMain.mat']);
codeSet = tempLoad.('pairsSoFar');
%% Simulation properties

% Define codes used for excitation
codes = cell(0);
%codesToUse = [3     7    13    15    53    55    63    95   119];%   125   167   181];
codesToUse = [1];% 7 9 11 13 15 17];

for i = 1:length(codesToUse)
    codes{i}.code = codeSet(codesToUse(i), :);
    codes{i}.ccode = codeSet(codesToUse(i)+1, :);
end

no_lines = 251;

% Transducer properties
f0 = 6.67e6;              %  Central frequency                        [Hz]
number_cycles = 2;
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     

width           = 0.2/1000;       % Width of element [m]
height  = 5/1000;         % Height of element [m]
kerf            = 0.02/1000;      % Kerf [m]

lambda = c / f0;       % Wavelength                                [m]
%pitch = lambda / 2;    % Pitch - center-to-center                  [m]
%width = .95*pitch;     % Width of the element                      [m]
%kerf = pitch - width;  % Inter-element spacing                     [m]
%height = 10/1000;      % Size in the Y direction                   [m]

rx_fnum         = 0;   % Receive f-number - Set to 0 to turn off
                       % constant F-number reconstruction
no_active = 192;        % Number of active elements on transmit and receive

%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:number_cycles/f0));
%impulse_response = impulse_response.*hanning(length(impulse_response))';
excitation = impulse_response;
%excitation = sin(2*pi*f0*(0:1/fs:number_cycles/f0));

%  Define the phantom
pht_pos = [%0 0 20;
           %0 0 30;
           0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
           -6 0 30; -6 0 35; -6 0 40; -6 0 45; -6 0 50;
           6 0 30; 6 0 35; 6 0 40; 6 0 45; 6 0 50;

           %-6 0 30; 0 0 30; 6 0 30;
           %-6 0 40; 0 0 40; 6 0 40;
           %3 0 20;
           %3 0 30;
           %6 0 40;
           %3 0 50;
           %-3 0 20;
           %-3 0 30;
           %-6 0 40;
           %-3 0 50;
           ] / 1000;         %  The position of the phantom
pht_pos = pht_pos;
%pht_pos(end+1, :) = [0 0 40]/1000;
pht_amp = ones(size(pht_pos, 1),1);      %  The amplitude of the back-scatter

% Calculate minimum and maximum samples and times for phantom
 [Rmax, Rmin, Tmin, Smin, max_code_length, Smin_c, Smax_c, no_rf_samples, no_rf_samples_c] =...
    calcSampleTimeRanges(pht_pos, codes, c, fs);

%  Initialize Field II and BFT
field_init(-1);
bft_init;

%  Set Field and BFT parameters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);
set_sampling(fs);
set_field('use_triangles', 0);

bft_no_lines(1);

% Create FIELD II transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Create BFT reconstruction array
xdc = bft_linear_array(no_active, width, kerf);

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

d_x_line = width*no_active / (no_lines-1);
x_line = -(no_lines-1) / 2 * d_x_line;
z_points = linspace(Rmin, Rmax, 100);
no_points = length(z_points);
T = z_points/c*2;
x_lines = x_line+d_x_line*(0:no_lines-1);

%% Simulate imaging with each code
% rf_data generated from imaging with each code and its associated focus
% zone is added together to simulate imaging with all codes. This is then 
% repeated for complementary versions.

% Rf data matrix holds non-complementary and complementary simulated
% RF data.
    
rf_data_m = zeros(no_rf_samples_c, no_elements, 2);
xmt_info = xdc_get(xmt, 'rect');
xElements = xmt_info(8, :);

no_active_tx = 64;
no_active_rx = 65;
no_focal_zones_x = 3;
focal_spacing_x = 6/1000;
bf_image = zeros(no_rf_samples, no_lines);

image_width = 28/1000;

focalZoneSpacing_x = ceil(no_lines/length(codes));
x = -image_width/2;
x_lines = linspace(-image_width/2, image_width/2, no_lines);
for lineNo = 1:focalZoneSpacing_x
    disp(['Transmit ', num2str(lineNo)]);
    % Build up the beam for each pair
    beam = zeros(1500, no_elements, 2);
    forceMaxDelay = 55;
    maxDelays = zeros(1, length(codes));
    for i = 1:length(codes)
        currentLine = lineNo+focalZoneSpacing_x*(i-1);
        if (currentLine > no_lines)
            continue;
        end
        
        focus = [x_lines(currentLine) 0 30/1000];
        % Find the apodization for this x focal zone
        N_pre_tx = round(focus(1)/(width+kerf) + no_elements/2 - no_active_tx/2);
        N_post_tx = no_elements - N_pre_tx - no_active_tx;
        apo_vector_tx = [zeros(1, N_pre_tx) ones(1, no_active_tx) zeros(1, N_post_tx)];
        
        N_pre_rx = round(focus(1)/(width+kerf) + no_elements/2 - no_active_rx/2);
        N_post_rx = no_elements - N_pre_rx - no_active_rx;
        apo_vector_rx = [zeros(1, N_pre_rx) ones(1, no_active_rx) zeros(1, N_post_rx)];
        
        tmp = find(apo_vector_tx);
        apoStart = tmp(1);
        apoEnd = tmp(end);

        codes{i}.lineNo = currentLine;
        codes{i}.apod_tx = apo_vector_tx;
        codes{i}.apod_rx = apo_vector_rx;
        
        focus(1) = 0;
        xEle = xElements(no_elements/2-no_active_tx/2+1:no_elements/2+no_active_tx/2);
        tmpBeam = focusBeam(excitation, focus, xEle, fs, c, forceMaxDelay);
        beam(1:size(tmpBeam, 1), apoStart:apoEnd, 1) = beam(1:size(tmpBeam, 1), apoStart:apoEnd, 1) + tmpBeam;

        %[tmpBeam, maxDelay] = focusBeam(codes{i}.ccode, focus, xEle, fs, c, forceMaxDelay);
        %beam(1:size(tmpBeam, 1), apoStart:apoEnd, 2) = beam(1:size(tmpBeam, 1), apoStart:apoEnd, 2) + tmpBeam;
    end
    beam = beam(1:forceMaxDelay+length(excitation), :, :);

    xdc_focus_times(xmt, 0, zeros(1, no_elements));
    xdc_focus_times(rcv, 0, zeros(1, no_elements));

    ele_waveform(xmt, (1:no_elements)', beam(:, :, 1)');
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 1) = rf_data;

%     ele_waveform(xmt, (1:no_elements)', beam(:, :, 2)');
%     [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
%     rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
%     rf_data_m(:, :, 2) = rf_data;

    rx_fnum = 0;
    %% Decode images
    % For each code, cross-correlate rf_data_sum with each code and add the
    % result to rf_data_decoded
    for i = 1:length(codes)
        % Beamform image for only lines we care about
        x_line = x_lines(codes{i}.lineNo);
        % Focus everywhere in the line
        bft_center_focus([x_line 0 0]);
        %bft_focus(xdc, 0, [x_line 0 40/1000]);
        bft_dynamic_focus(xdc, 0, 0); 

        % Constant F-Num reconstruction
        if (rx_fnum == 0)
            % Use all elements if fnum = 0

            bft_apodization(xdc, 0, ones(1, no_elements));
            %bft_apodization(xdc, 0, codes{i}.apod_rx);
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

            bft_apodization(xdc, T, apo_vector_rx);
        end
        bf_temp = bft_beamform(Tmin-forceMaxDelay/fs, rf_data_m(1:no_rf_samples, :, 1));
        bf_image(:, codes{i}.lineNo) = bf_temp;
    end
end
%bf_temp = bft_beamform(Tmin-maxDelays(i), sum(rf_data_decoded, 3));
%bf_image = bf_temp;
% Can either add decoded data for non-complementary and complementary then
% beamform.
%rf_data_decoded = rf_data_nc_decoded + rf_data_c_decoded;
%bf_temp = bft_beamform(Tmin, sum(rf_data_decoded, 3));

% Or we can beamform decoded data then add results
%bf_temp = bft_beamform(Tmin, rf_data_nc_decoded);
%bf_temp = bf_temp + bft_beamform(Tmin, rf_data_c_decoded);

field_end
bft_end

% Perform envelope detection and normalize
env_bf = abs(hilbert(bf_image));
env_bf = env_bf / max(max(env_bf));

%% Plot image
figure;
imagesc(x_lines*1000, ([Rmin Rmax]-forceMaxDelay/fs*c/2)*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);
