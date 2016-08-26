% compCodesWalkingMultiple.m

% This program simulates point spread functions due to pulse-echo B-mode 
% imaging using a Siemens Antares ultrasound system.

% By Roger Zemp. Updated 9 Mar 2011.
% Form three lines at once, each using a different pair of complementary
% codes. Then, form the next three lines on the next pair of transmits.

% Better point spread function than other methods we used
% Some stiching issues, some skew issues.

%% Simulation properties

% Load codes used for excitation
tempLoad = load(['../../Complementary Pairs/compPairs_len_10_simMain.mat']);
codeSet = tempLoad.('pairsSoFar');

% Store codes in a cell array for later use
codes = cell(0);
codesToUse = [1]; % Which codes to use from code set
for i = 1:length(codesToUse)
    codes{i}.code = codeSet(codesToUse(i), :); % 1st pair code
    codes{i}.ccode = codeSet(codesToUse(i)+1, :); % 2nd pair code
    codes{i}.focusZ = 30/1000;
end

% Transducer properties
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
fs = 100e6;             % Sampling frequency                       [Hz]
c = 1540;               % Speed of sound                           [m/s]
no_elements = 192;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 10/1000;       % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 

no_lines = 251;         % Number of lines in image
no_active_tx = 64;      % Number of active elements for a transmit 
                        % sub-aperture
rx_fnum = 2.1;          % Receive f-number for calculating number of 
                        % active elements for receive sub-aperture
rx_focus = [0 0 0.030]; % Receive focus for calculaing number of active
                        % elements for receive sub-aperture

image_width = 28/1000;  % Image width to simulate [m]

forceMaxDelay = 55;

%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:number_cycles/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
pht_pos = [           
           0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
           -6 0 30; -6 0 35; -6 0 40; -6 0 45; -6 0 50;
           6 0 30; 6 0 35; 6 0 40; 6 0 45; 6 0 50;
           ] / 1000; %  The position of the phantom
pht_amp = ones(size(pht_pos, 1),1); %  The amplitude of the back-scatter

% Starting and ending distance to image
[Rmax, Rmin, Tmin, Smin, max_code_length, ...
 Smin_c,Smax_c, no_rf_samples, no_rf_samples_c] = ...
calcSampleTimeRanges(pht_pos, codes, c, fs);

%  Initialize Field II and BFT
field_init(-1);
bft_init;

%  Set Field and BFT parameters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);

bft_no_lines(1); % Number of lines to beamform at a time

% Create FIELD II transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Create BFT reconstruction array
xdc = bft_linear_array(no_elements, width, kerf);

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

% Set the apodization
xdc_apodization(xmt, 0, ones(1,no_elements))
xdc_apodization(rcv, 0, ones(1,no_elements))

% Set center for transmit and receive focusing
xdc_center_focus(xmt,[0 0 0])
xdc_center_focus(rcv,[0 0 0])

% Turn off Field II receive focusing
xdc_focus_times(rcv, 0, zeros(1,no_elements)); 

% Calculate number of active receive elements for a sub-aperture
rx_ap = rx_focus(3)/rx_fnum;
no_active_rx = round(rx_ap/(width+kerf));   

%% Simulate imaging
% For each set of lines, we image with each code and its complementary pair

% Store raw RF data for each element in rf_data_m
% 3rd dimension is transmit event index
rf_data_m = zeros(no_rf_samples_c, no_elements, 2);

% Get x position of transducer elements directly from Field II
xmt_info = xdc_get(xmt, 'rect');
xElements = xmt_info(8, :);

% Calculate focal zone spacing in x in units of lines
focalZoneSpacing_x = ceil(no_lines/length(codes));

% Calculate x positions of lines
x_lines = linspace(-image_width/2, image_width/2, no_lines);

% Stores beamformed image
bf_image = zeros(no_rf_samples, no_lines);
for lineNo = 1:focalZoneSpacing_x
    disp(['Line set ', num2str(lineNo)]);
    
    % Build up the beam for each pair of codes and sum the beams into one
    % 3rd dimension is transmit index for this line    
    beam = zeros(forceMaxDelay+max_code_length, no_elements, 2);
    
    % Go through each set of code pairs
    for i = 1:length(codes)
        % Get index of current line
        currentLine = lineNo + focalZoneSpacing_x*(i-1);
        if (currentLine > no_lines) % Stop once beamformed all lines
            continue;
        end
        
        % Set transmit focus for current line
        focus = [x_lines(currentLine) 0 codes{i}.focusZ];
        
        % Calculate transmit apodization        
        N_pre_tx = round(focus(1)/(width+kerf) + no_elements/2 - no_active_tx/2);
        N_post_tx = no_elements - N_pre_tx - no_active_tx;
        apo_vector_tx = [zeros(1, N_pre_tx) ones(1, no_active_tx) zeros(1, N_post_tx)];
        
        % Calculate receive apodization
        N_pre_rx = round(focus(1)/(width+kerf) + no_elements/2 - no_active_rx/2);
        N_post_rx = no_elements - N_pre_rx - no_active_rx;
        apo_vector_rx = [zeros(1, N_pre_rx) ones(1, no_active_rx) zeros(1, N_post_rx)];
        
        % Storing current line number and apodization
        codes{i}.lineNo = currentLine;
        codes{i}.apod_tx = apo_vector_tx;
        codes{i}.apod_rx = apo_vector_rx;
        
        % Find first and last elements in transmit apodization
        tmp = find(apo_vector_tx);
        apoStart = tmp(1);
        apoEnd = tmp(end);
        
        % Manually transmit focus @focus for each code in the pair
        xEle = xElements(apoStart:apoEnd); % Get x position of elements 
                                           % in transmit apodization        
        tmpBeam = focusBeam(codes{i}.code, focus, xEle, fs, c, forceMaxDelay);
        beam(:, apoStart:apoEnd, 1) = beam(:, apoStart:apoEnd, 1) + tmpBeam;

        tmpBeam = focusBeam(codes{i}.ccode, focus, xEle, fs, c, forceMaxDelay);
        beam(:, apoStart:apoEnd, 2) = beam(:, apoStart:apoEnd, 2) + tmpBeam;
    end
    
    % Don't tranmit or receive focus using Field II
    xdc_focus_times(xmt, 0, zeros(1, no_elements));
    xdc_focus_times(rcv, 0, zeros(1, no_elements));
    
    % Set transducer excitation and simulate with calc_scat_multi
    % And align scatter data (force data to lie in a specified time window)
    % The corresponding sample window is [Smin_c, Smax_c]
    ele_waveform(xmt, (1:no_elements)', beam(:, :, 1)'); 
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 1) = rf_data;
    
    % Repeat for second code in pair
    ele_waveform(xmt, (1:no_elements)', beam(:, :, 2)');
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 2) = rf_data;

    %% Decode each line in the current line set for this transmit event

    % Store RF data for a single code pair (3rd dimension is transmit)
    rf_data_decoded = zeros(no_rf_samples, no_elements, 2);
    for i = 1:length(codes) % Loop through each code pair      
        
        % Decode with respect to first code in current pair
        temp_decoded = conv2(rf_data_m(:, :, 1), rot90(conj(codes{i}.code'), 2), 'valid');
        rf_data_decoded(:, :, 1) = temp_decoded(1:no_rf_samples, :);
        
        % Decode with respect to second code in current pair
        temp_decoded = conv2(rf_data_m(:, :, 2), rot90(conj(codes{i}.ccode'), 2), 'valid');
        rf_data_decoded(:, :, 2) = temp_decoded(1:no_rf_samples, :);
        
        % Beamform image for the current line
        
        % Set center focus for current line and enable dynamic focusing
        x_line = x_lines(codes{i}.lineNo);
        bft_center_focus([x_line 0 0]);       
        bft_dynamic_focus(xdc, 0, 0); 

        bft_apodization(xdc, 0, ones(1, no_elements));
        %bft_apodization(xdc, 0, codes{i}.apod_rx);
        bf_temp = bft_beamform(Tmin-forceMaxDelay/fs, sum(rf_data_decoded, 3));
        bf_image(:, codes{i}.lineNo) = bf_temp;
    end
end

% Release memory used by Field II and BFT
field_end
bft_end

% Perform envelope detection and normalize
env_bf = abs(hilbert(bf_image));
env_bf = env_bf / max(max(env_bf));

%% Plot image
figure;
imagesc(x_lines*1000, ([Rmin Rmax]-(forceMaxDelay/fs*c/2))*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);
