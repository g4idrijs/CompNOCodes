% walkAp_compl_sim
% Simulate a flash image with complementary codes at many focal zones
% VERSION 1.0, August 10, 2016 (Tarek)
%
% This version uses BFT to beamform images
% 
% Version 1.1, August 11, 2016 (David)
% Adds calcuation of axial and lateral full width half maximum
% as well as signal to noise ratio
%
% Currently no support for:
% -Constant f number
% -Signal to noise ratio
% Currently working on looping over aperatures instead of codes

%% Initial setup (ex. transducer, phantom, image range properties) 
clear

% Add white Gaussian noise
useNoise = 0;

% Define codes used for excitation
% along with point for transmit focusing (center of focus domain)
codes = cell(0);
codes{1}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
codes{1}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
codes{1}.focus = [0/1000 0 50/1000]; 

% codes{2}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
% codes{2}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
% codes{2}.focus = [0/1000 0 50/1000]; 
% 
% codes{3}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
% codes{3}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
% codes{3}.focus = [10/1000 0 50/1000]; 

% codes{3}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
% codes{3}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
% codes{3}.focus = [10/1000 0 50/1000]; 



% tempDat = load('bestPairsSQP41.mat');
% x = tempDat.('x');
% codes{1}.code = x(1,:); 
% codes{1}.ccode = x(2,:);

% Transducer properties
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     
no_active_tx = 64;      %  Number of transmitting elements in walking aperture
no_active_rx = 64;      %  Number of receiving elements in walking apertur

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]

% rx_fnum         = 0; % 2.1;   % Receive f-number - Set to 0 to turn off
                       % constants F-number reconstruction

%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:1/f0)); 
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
pht_pos = [0 0 50/1000]; % Position
pht_amp = 20; % Strength of scattering

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

%  Initialize Field II and beamforming toolbox (BFT)
field_init(-1);
bft_init;

%  Set Field and BFT parameters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);

% Create transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Create beamforming array
xdc = bft_linear_array(no_elements, width, kerf);

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

%% Imaging preparation (ex. line setup)

% Don't receive focus (will do that using beam forming toolbox)
xdc_focus_times(rcv, 0, zeros(1,no_elements));

% One focus region per complementary pair
no_CodePairs = length(codes);
no_FocusDomains = no_CodePairs; 

% Each focus domain has a number of lines associated with it
no_linesPerDomain = 100;
no_lines = no_linesPerDomain*no_FocusDomains; % Total number of lines

% Calculate line widths assuming constant focus domain spacing
% This also assumes the codes are listed in order of focus region
widthFocusDomain = no_active_rx*width; %10/1000; %abs(codes{2}.focus(1)  - codes{1}.focus(1));
d_x_line = widthFocusDomain/no_linesPerDomain;

 %% Imaging loop

% We have to process each transmit separately
firstCodeCorr = zeros(no_rf_samples, no_elements,no_FocusDomains); % Correlate with codes on first transmit
secCodeCorr = zeros(no_rf_samples, no_elements,no_FocusDomains);  % Correlate with codes on second transmit
for currTransmit = 1:2
    % Get RF data from each focus domain   
    % Each column has data from one transducer element
    RFTransRunTot = zeros(no_rf_samples_c,no_elements);  
    for currFocDomain = 1:no_FocusDomains
        % Set transmit focus at center of focus domain
        transFocus = codes{currFocDomain}.focus;
        xFocus = transFocus(1);     

        % Set apodization centered around center of focus domain
        % (for transmit - we are forced to receive on all since all
        % code pairs are fired at the same time)
        [apo_vector_tx,apo_vector_rx] = getApodization(xFocus,no_active_tx,no_active_rx,width, kerf, no_elements);
        xdc_apodization(xmt, 0, apo_vector_tx);
        xdc_apodization(rcv, 0, ones(1,no_elements));

        % Get the current complementary code pair
        currCodes = codes{currFocDomain};

        % Set transmit focus at center of focus domain
        xdc_center_focus(xmt,[xFocus, 0, 0]) % Still not clear on this line..
        xdc_focus(xmt, 0, currCodes.focus);

        % Get RF scattering data by firing...    
        % ...first code in pair
        if(currTransmit == 1)
            xdc_excitation(xmt, currCodes.code);
            [codeRF, start_timeCode] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp); 
        else
            %...second code in pair
            xdc_excitation(xmt, currCodes.ccode);
            [codeRF, start_timeCode] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp); 
        end
        
        % Align received RF data in time
        codeRF = alignRF( codeRF, start_timeCode,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);           

        % Add new RF data to already received RF data (for this transmit)
        RFTransRunTot = RFTransRunTot + codeRF;
    end
     
    % We have collected the data for one transmit
    % Decode to get desired data for each focus domain
    for currFocDomain = 1:no_FocusDomains
        % Get current code to cross correlate with for this focus domain
        if (currTransmit == 1)
            currFocDomCode = codes{currFocDomain}.code;
        else
            currFocDomCode = codes{currFocDomain}.ccode;   
        end 
        % Cross correlate
        currFocDomDec = conv2(RFTransRunTot, rot90(conj(currFocDomCode'), 2), 'valid');
        currFocDomDec = currFocDomDec(1:no_rf_samples, :);
        
        % Store results
        if (currTransmit == 1)
            firstCodeCorr(:,:,currFocDomain) = currFocDomDec;
        else
            secCodeCorr(:,:,currFocDomain) = currFocDomDec;            
        end
    end      
end
 
% Sum the results of cross correlating with the different codes
% Each matrix in 3D dimension is the decoded data for a focus domain
% (has signal received by each transducer - the signal is that relevant
% to that focus domain)
sumDec =  firstCodeCorr + secCodeCorr;

% Beamform each focus domain
bfImag = zeros(no_rf_samples,no_lines); % Just store one line per focus domain for now
lineCount = 1;
for currFocDomain = 1:no_FocusDomains
    % Pick the center receive focus location for this focus domain
    % Beamform along the several lines in the focus domain
    centFocus = codes{currFocDomain}.focus;
    xCentFocus = centFocus(1);
    
    % Set focusing for beamforming
    
    % Set beamforming apodization
    [apo_vector_tx,apo_vector_rx] = getApodization(xFocus,no_active_tx,no_active_rx,width, kerf, no_elements);
    bft_apodization(xdc, 0, apo_vector_rx); 
    
    % Starting x position of receive focus
    currX = xCentFocus - widthFocusDomain/2;
    
    for currLine = 1:no_linesPerDomain
        % Set current receive focus        
        bft_center_focus([currX 0 0]); % Set up start of line for dynamic focus       
        bft_dynamic_focus(xdc, 0, 0);  % Set direction of dynamic focus
        %bft_focus(xdc, 0, [currX 0 0])
        
        % Beamform and store resulting line
        currLineBf = bft_beamform(Tmin, sumDec(:,:,currFocDomain));     
        bfImag(:,lineCount) = currLineBf;
        
        % Move receive focus
        currX = currX + d_x_line;
        
        % Keep track which line we are on
        lineCount = lineCount + 1;
    end
end

% Normalize and do envelope detection
env_bf = abs(hilbert(bfImag));
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
 
%% ALL BELOW HERE IS OLD
%% We transmit the first code in each pair, then the second in each pair
no_Transmits = 2;
transmitResults = cell(2,1);
transmitResults{1} = zeros(no_rf_samples,no_lines);
transmitResults{2} = zeros(no_rf_samples,no_lines);

for currTransmit = 1:no_Transmits
    % For each transmit, we image and beamform each line
    for currLine = 1:no_lines
        % Calculate and set apodization for this line
        xFocus = x_start + (currLine-1)*d_x_line; % Focus at [xFocus, 0, 0] for this line
        [apo_vector_tx,apo_vector_rx] = getApodization(xFocus,no_active_tx,no_active_rx,width, kerf, no_elements);
        xdc_apodization(xmt, 0, apo_vector_tx);
        xdc_apodization(rcv, 0, apo_vector_rx);
        
        % Set center reference point for transmit focus for this line
        xdc_center_focus(xmt,[xFocus, 0, 0])
        
        % For each line, we fire each relevant and store the results
        % Currently all codes in this transmit are fired for each line       
        RF_runningTot = 0; % Add the results from each code to this total
        for currCode = 1:no_CodePairs  
            % Get current code, and set as transducer excitation
            if (currTransmit == 1)
                code = codes{currCode}.code;
            else
                code = codes{currCode}.ccode;   
            end            
            xdc_excitation(xmt, code);
            
            % Set the transmit focus for this code
            xdc_focus(xmt, 0, codes{currCode}.focus);
            
            % Get RF data for transducers that fire
            [currCodeRF, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp); 
            
             % Make sure the RF data exists exactly during the desired times
            currCodeRF = alignRF(currCodeRF,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);           
                        
            % Add to the RF running total for this line
            RF_runningTot = RF_runningTot + currCodeRF;
        end
        
        % Decode the data for this line
        % (cross correlate with each code in this transmit and add to a running total)
        currLineDec = 0;
        for currCode = 1:no_CodePairs
            % Get current code
            if (currTransmit == 1)
                code = codes{currCode}.code;
            else
                code = codes{currCode}.ccode;   
            end 
            % Cross correlate with current code
            currCodeDec = conv2(RF_runningTot, rot90(conj(code'), 2), 'valid');
            currCodeDec = currCodeDec(1:no_rf_samples, :);  
            
            % Add results to running total (over all codes used)
            currLineDec = currLineDec + currCodeDec;
        end
        
        % Beamform this line
        bft_center_focus([xFocus 0 0]); % Set up start of line for dynamic focus       
        bft_dynamic_focus(xdc, 0, 0);  % Set direction of dynamic focus
        bft_apodization(xdc, 0, apo_vector_rx); % Beamforming apodization
        currLineBf = bft_beamform(Tmin, currLineDec); % Beamform
        
        % Add the beamformed line to the results for this transmit
        transmitResults{currTransmit}(:,currLine) = currLineBf + transmitResults{currTransmit}(:,currLine);
    end
end

% Add the transmit matrices to create a beamformed imaged
% (not normalized, not enveloped)
bfImag = plus(transmitResults{:});

% Normalize and do envelope detection
env_bf = abs(hilbert(bfImag));
env_bf = env_bf / max(max(env_bf));

%% Plot image
figure;
imagesc([-1/2 1/2]*no_lines*d_x_line*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);

%% Set up beamforming
for i = 1 : no_lines
   % Set reference point for dynamic focusing for this line
   % Dynamic focusing focuses everywhere in this line
    bft_center_focus([x_line 0 0], i); 
    
    % Set angle of dynamic focus line for this line
    % We set it to point straight at the focus
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
for i = 1:2*length(codes) % Go through both first code and second code for each pair
        
    % Extract code and focus information for current code
    if (i > length(codes)) % Working with second code "ccode"
        code = codes{i-length(codes)}.ccode;
        focus = codes{i-length(codes)}.focus;
        col = 2; % Where to store results
    else
        code = codes{i}.code; % Working with first code "code"
        focus = codes{i}.focus;
        col = 1; % Where to store results
    end        
    
    % Set focus for this particular code
    % (each code has its own focus)
    xdc_focus(xmt, 0, focus);
    xdc_excitation(xmt, code);
    
    % Get RF_data from each transducer that fired and align it in time
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

% Add noise
if(useNoise == 1)
    % Add white gaussian noise (prior to beamforming)     
    noisePowerdBW =  -4.45e2; % Use -inf to turn off noise
    rf_data_m(:,:,1)  = rf_data_m(:,:,1) + wgn(size(rf_data_m,1), size(rf_data_m,2),noisePowerdBW);
    rf_data_m(:,:,2)  = rf_data_m(:,:,2) + wgn(size(rf_data_m,1), size(rf_data_m,2),noisePowerdBW);
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

%% Calculate signal to noise ratio
if(useNoise == 1)
    % Signal / sd(noise)
    % std(noise) =  about 6e-21 with Gaussian white noise at -4.45e2 dBW.
    % and with the following codes:
    %     codes{1}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
    %     codes{1}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
    % (see the findNoise script)
    maxSig = max(max(bf_temp));
    stdNoise = 6e-21;
    sigToNoise = maxSig/stdNoise;
    disp('Signal to noise ratio:');
    disp(sigToNoise);
end

%% Perform envelope detection and normalize
env_bf = abs(hilbert(bf_temp));
env_bf = env_bf / max(max(env_bf));

%% Calculate full width half maximum resolution
% (Needs generalization to work with a several point phantom)
% Find the index of the maximum signal
[maxVal,maxLocLin] = max(env_bf(:));
[maxLocZ, maxLocX] = ind2sub(size(env_bf),maxLocLin);

% Axial (z direction)
zIntVect = env_bf(:,maxLocX);
fwhmIndicesZ = fullWidthHM(zIntVect); % FWHM in indices
lenPerIndexZ = (Rmax - Rmin)/numel(zIntVect);
fwhmLengthZ = lenPerIndexZ*fwhmIndicesZ; % Convert to length
disp('Axial full width half maximum: (mm)')
disp(fwhmLengthZ*1000)

% Lateral (x direction)
xIntVect = env_bf(maxLocZ,:);
fwhmIndicesX = fullWidthHM(xIntVect); % FWHM in indices
lenPerIndexX = (no_lines*d_x_line)/numel(xIntVect);
fwhmLengthX = lenPerIndexX*fwhmIndicesX; % Convert to length
disp('Lateral full width half maximum: (mm)')
disp(fwhmLengthX*1000)

%% Plot image
figure;
imagesc([-1/2 1/2]*no_lines*d_x_line*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);