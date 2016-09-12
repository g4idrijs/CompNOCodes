% Copy of compCodesWalkingMultiple (September 8, 2016)
% Purpose: try out calc_hp: visualize transducer output

% compCodesWalkingMultiple.m

% Authors: Tarek, David

% Description Uses Field II and BFT to simulate a walking aperture of 
% parallel focal zones in both X and Z directions using complementary 
% codes. Number of parallel zones in X direction determine how many lines 
% will be simultaneously walked, and number of parallel zones in Z 
% direction determine how many parallel focal zones in the Z direction.

clear
addpath('../../Field_II', '../../bft_64bit');

%% Simulation properties

% Load codes used for excitation
% tempLoad = load(['../../../Complementary Pairs/compPairs_len_10_simMain.mat']);
% codeSet = tempLoad.('pairsSoFar');
tempLoad =  load('../../../Complementary Pairs/len10_100codes_minInt2.mat');
codeSet = tempLoad.('x');

% tempLoad =  load('../../../Complementary Pairs/len5_3codes.mat');
% codeSet = tempLoad.('x');

% tempLoad1 = load(['../../../Complementary Pairs/len10_4codes_set1.mat']);
% tempLoad2 = load(['../../../Complementary Pairs/len10_4codes_set2.mat']);
% tempLoad3 = load(['../../../Complementary Pairs/len10_4codes_set3.mat']);
% tempLoad4 = load(['../../../Complementary Pairs/len10_4codes_set4.mat']);
% tempLoad = [tempLoad1.('x'); tempLoad2.('x'); tempLoad3.('x'); tempLoad4.('x')];
% codeSet = tempLoad;

% Store codes in a cell array for later use
numCodesX = 1; % Number of parallel focal zones in X
numCodesZ = 1; % Number of parallel focal zones in Z

codes = cell(0);
codesToUse = [1:2:numCodesX*numCodesZ*2-1]; % Which codes to use from code set

if(numCodesZ == 1)
   focalPoints_z = 40/1000; 
else
   focalPoints_z = linspace(30, 50, numCodesZ)/1000;
end
 
for i = 1:numCodesX
    for j = 1:numCodesZ
        codes{i}.code{j} = [1]; %codeSet(codesToUse(1+(i-1)*numCodesZ+(j-1)), :); % 1st pair code
        codes{i}.ccode{j} = [1]; %codeSet(codesToUse(1+(i-1)*numCodesZ+(j-1))+1, :); % 2nd pair code
        codes{i}.focusZ(j) = focalPoints_z(j);
    end
end

% Transducer properties
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
fs = 100e6;             % Sampling frequency                       [Hz]
c = 1540;               % Speed of sound                           [m/s]
no_elements = 192;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 

no_lines = 151;         % Number of lines in image 
no_active_tx = 64;     % Number of active elements for transmit 
                        % sub-aperture
rx_fnum_constant = 0;   % F-number for constant f-num reconstruction.
                        % If this is set to 0 then no constant f-num
                        % reconstrution is done and rx_fnum will be used
                        % instead to determine number of active elements
                        % to use on receive.
rx_fnum = 0;          % Receive f-number for calculating number of 
                        % active elements for receive sub-aperture. Value
                        % of 0 means all elements will be used for
                        % reconstruction. Need to also set rx_focus to 
                        % calculate appropriate aperture.
rx_focus = [0 0 40/1000]; % Receive focus for calculaing number of active
                        % elements for receive sub-aperture [m]
                        
rxApodFunc = @(x) hanning(x); % Function to use for receive apodization

image_width = 12/1000;  % Image width to simulate [m]

forceMaxDelay = 50;   % Forces the maximum delay in manual focusing to 
                        % this value. This is needed to align multiple 
                        % focused beams in the Z direction. [samples]

%  Define the impulse response of the transducer
Ts = 1 /fs; % Sampling period
impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%% Extend codes
% Insert zeroes in codes to stop impulse responses from interfering
% numZeroInsert = numel(impulse_response) - 1;
% for i = 1:numCodesX*numCodesZ
%     N = numel(codes{i}.code{1});
%     
%     codes{i}.code{1} = matintrlv([codes{i}.code{1}  zeros(1, N*numZeroInsert);],numZeroInsert+1,N);
%     codes{i}.ccode{1} = matintrlv([codes{i}.ccode{1}  zeros(1, N*numZeroInsert);],numZeroInsert+1,N);
% end

% Extend codes by repeating elements
timesToRepeat = 20; %numel(impulse_response) - 1;
for i = 1:numCodesX*numCodesZ    
    
    codes{i}.code{1} = repelem(codes{i}.code{1} ,1,timesToRepeat);
    codes{i}.ccode{1} = repelem(codes{i}.ccode{1} ,1,timesToRepeat);
end

%% Phantom definition

%  Define the phantom
pht_pos = [     0 0 40  
%            0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
%            -6 0 30; -6 0 35; -6 0 40; -6 0 45; -6 0 50;
%            6 0 30; 6 0 35; 6 0 40; 6 0 45; 6 0 50;
           ] / 1000; %  The position of the phantom
pht_amp = ones(size(pht_pos, 1),1); %  The amplitude of the back-scatter

% Calculate start and end distance, sample, and time for phantom
[Rmax, Rmin, Tmin, Smin, max_code_length, ...
 Smin_c,Smax_c, no_rf_samples, no_rf_samples_c] = ...
calcSampleTimeRanges(pht_pos, codes, c, fs);

%% Field II and BFT setup
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

% Define set of z points and time points that will be used for constant
% f-number reconstruction
z_points = linspace(Rmin, Rmax, 100);
no_points = length(z_points);
T = z_points/c*2;

%% Simulate imaging
% For each set of lines, we image with each code and its complementary pair

% Store raw RF data for each element in rf_data_m
% 3rd dimension is transmit event/code pair index
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
outputMsg = '';
for lineNo = 1:focalZoneSpacing_x
    % Print current line number on same line
    fprintf(repmat('\b', 1, length(outputMsg)));
    outputMsg = sprintf('Line set %d/%d', lineNo, focalZoneSpacing_x);
    fprintf(outputMsg);
    
    % Build up the beam for each pair of codes and sum the beams into one
    % 3rd dimension is transmit index for this line    
    beam = zeros(forceMaxDelay+max_code_length, no_elements, 2);
    
    % Go through each set of code pairs
    for i = 1:length(codes)
        % Get index of current line
        currentLine = lineNo + focalZoneSpacing_x*(i-1);
        if (currentLine > no_lines) % Stop once beamformed all lines
            codes{i}.lineNo = 0;
            continue;
        end
        xLine = x_lines(currentLine);
        
        % Calculate transmit apodization from no_active_tx
        if no_active_tx == no_elements
            apo_vector_tx = ones(1, no_elements);
        else
            no_active_tx_tmp = no_active_tx;
            N_pre_tx = round(xLine/(width+kerf) + no_elements/2 - no_active_tx/2);
            if (N_pre_tx < 0)
                no_active_tx_tmp = no_active_tx_tmp + N_pre_tx;
                N_pre_tx = 0;
            end
            N_post_tx = no_elements - N_pre_tx - no_active_tx_tmp;
            if (N_post_tx < 0)
                no_active_tx_tmp = no_active_tx_tmp + N_post_tx;
                N_post_tx = 0;
            end
            apo_vector_tx = [zeros(1, N_pre_tx) ones(1, no_active_tx_tmp) zeros(1, N_post_tx)];
        end
        
        % Calculate receieve apodization from rx_fnum
        if ~rx_fnum
            % If rx_fnum = 0 then use all elements for receive apodization
            apo_vector_rx = ones(1, no_elements);
        else
            % Calculate receive apodization for rx_fnum and rx_focus 
            N_pre_rx = round(xLine/(width+kerf) + no_elements/2 - no_active_rx/2);
            N_post_rx = no_elements - N_pre_rx - no_active_rx;
            apo_vector_rx = [zeros(1, N_pre_rx) ones(1, no_active_rx) zeros(1, N_post_rx)];
        end
        % Store current line number and apodization into current code
        codes{i}.lineNo = currentLine;
        line(i, lineNo) = currentLine;
        
        xxline(i, lineNo) = xLine;
        codes{i}.apod_tx = apo_vector_tx;
        codes{i}.apod_rx = apo_vector_rx;
        
        % Find first and last elements in transmit apodization
        tmp = find(apo_vector_tx);
        apoStart = tmp(1);
        apoEnd = tmp(end);

        % Manually transmit focus @focus for each code in the pair
        xEle = xElements(apoStart:apoEnd); % Get x position of elements 
                                           % in transmit apodization
            
        for j = 1:length(codes{i}.code)
            % Set transmit focus for current line
            focus = [xLine 0 codes{i}.focusZ(j)];
        
            tmpBeam = focusBeam(codes{i}.code{j}, focus, xEle, fs, c, forceMaxDelay);
            beam(:, apoStart:apoEnd, 1) = beam(:, apoStart:apoEnd, 1) + tmpBeam;

            tmpBeam = focusBeam(codes{i}.ccode{j}, focus, xEle, fs, c, forceMaxDelay);
            beam(:, apoStart:apoEnd, 2) = beam(:, apoStart:apoEnd, 2) + tmpBeam;
        end
    end
    
    % Don't tranmit or receive focus using Field II
    xdc_focus_times(xmt, 0, zeros(1, no_elements));
    xdc_focus_times(rcv, 0, zeros(1, no_elements));
    
    % Set transducer excitation and simulate with calc_scat_multi
    % Align scatter data (force data to lie in a specified time window)
    % The corresponding sample window is [Smin_c, Smax_c]
    % Also need to shift rf_data by forceMaxDelay before or after aligning 
    % to account for focusing delays.
    ele_waveform(xmt, (1:no_elements)', beam(:, :, 1)');     
   
    %% Try visualizing beam
    if(xLine == 0)
        % Plot excitation over time at a single points
        plotSingLoc = 0;
        if(plotSingLoc == 1)
            visLoc = [0 0 60];
            [h,start_time] = calc_hp(xmt,visLoc/1000);
            figure
            plot(h)
            title(['Energy at [', num2str(visLoc),'] mm']);
        end

        % Get excitation at a range of positions over time
        xRange = [-10:0.1:10]/1000; % m
        zRange = [0:0.1:40]/1000; %[35:0.1:45]/1000; % m    
        pointsInterest = zeros(numel(zRange),3);
        pointsInterest(:,3) = zRange'; 

        % Determine indices of samples to keep
        [Rmax, Rmin, ~, ~, ~, Smin_cVis, Smax_cVis, ~, no_rf_samples_cVis] =...
        calcSampleTimeRanges(pointsInterest, codes, c, fs);
        Smin_cVis = round(Smin_cVis/2); % Adjust for no echo here
        Smax_cVis = round(Smax_cVis/2);
        no_rf_samples_cVis = round(no_rf_samples_cVis /2);

        for xRangeInd = 1:numel(xRange)
           disp(xRangeInd);
           for zRangeInd = 1:numel(zRange)
               currX = xRange(xRangeInd);
               currZ = zRange(zRangeInd);
               [currResp,start_timeVis] = calc_hp(xmt,[currX 0 currZ]);  
               
               currResp = alignRF(currResp,start_timeVis,fs,Smin_cVis,Smax_cVis,no_rf_samples_cVis,1);

               % Store sequence at right (x,z) location
               respMat(zRangeInd,xRangeInd,:) = currResp;
           end
        end

        %plot(currResp)   

        %% Show the animation
        figure
        numFrames = size(respMat,3);
        for frame = 300 %1:5:numFrames 
           imagesc([min(xRange), max(xRange)]*1000,[Rmin, Rmax]*1000, respMat(:,:,frame));

           title(['First Code Set. Line: ', num2str(lineNo),'. Frame: ', num2str(frame), ' of ', num2str(numFrames)]);
           colorbar;
           caxis([-0.1 0.1]*1e-11)
           ylabel('z (mm)')
           xlabel('x (mm)')

           pause(0.1);
        end

        %pause  
    end
    
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = [rf_data(1+forceMaxDelay:end, :); zeros(forceMaxDelay, size(rf_data, 2))];
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 1) = rf_data;
    
    % Repeat for second code in pair
    ele_waveform(xmt, (1:no_elements)', beam(:, :, 2)');    

    
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = [rf_data(1+forceMaxDelay:end, :); zeros(forceMaxDelay, size(rf_data, 2))];
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 2) = rf_data;    

    
    %% Decode each line in the current line set for this transmit event

    % Store RF data for a single code pair (3rd dimension is transmit)
    rf_data_decoded = zeros(no_rf_samples, no_elements, 2);
    for i = 1:length(codes) % Loop through each code pair
        if (codes{i}.lineNo < 1)
            continue;
        end
        % Set center focus for current line and enable dynamic focusing
        x_line = x_lines(codes{i}.lineNo);
        bft_center_focus([x_line 0 0]);       
        bft_dynamic_focus(xdc, 0, 0); 
        
        % Constant F-Num reconstruction
        if (~rx_fnum_constant)
            % Use apodization calculated during beam construction
            tmp = find(codes{i}.apod_rx);
            apodStart = tmp(1);
            apodEnd = tmp(end);
            numEndZero = length(codes{i}.apod_rx) - apodEnd;
            apod = [zeros(1,apodStart-1) rxApodFunc(apodEnd - apodStart + 1)' zeros(1,numEndZero)];
            
            bft_apodization(xdc, 0, apod);
        else
            % Use constant f-num reconstruction
            rx_ap_c           = z_points/rx_fnum_constant;
            no_active_rx_c    = min(round(rx_ap_c/(width+kerf)), no_elements);

            mid_elements = round(x_line/(width+kerf)+no_elements/2);
            start_elements = round(max(mid_element - no_active_rx_c/2 + 1, 1));
            end_elements = round(min(mid_element+no_active_rx_c/2, no_elements));

            apo_vector_rx = zeros(no_points, no_elements);
            for j = 1:no_points
                apo_vector_rx(j, start_elements(j):end_elements(j)) = rxApodFunc(end_elements(j)-start_elements(j)+1);
            end

            bft_apodization(xdc, T, apo_vector_rx);
        end
        
        for j = 1:length(codes{i}.code)
%             % Try correlating to codes after have gone through transducer
%             % (DEBUG)
%             scrCode = rand(1,2)*2-1; %conv(conv(codes{i}.code{j}, impulse_response),impulse_response);
%             scrCCode = rand(1,2)*2-1; %conv(conv(codes{i}.ccode{j}, impulse_response),impulse_response);
% 
%              
%             temp_decoded = conv2(rf_data_m(:, :, 1), rot90(conj(scrCode'), 2), 'valid');
%             rf_data_decoded(:, :, 1) = temp_decoded(1:no_rf_samples, :);
% 
%             temp_decoded = conv2(rf_data_m(:, :, 2), rot90(conj(scrCCode'), 2), 'valid');
%             rf_data_decoded(:, :, 2) = temp_decoded(1:no_rf_samples, :);
            
            % Decode with respect to first code in current pair   
            temp_decoded = conv2(rf_data_m(:, :, 1), rot90(conj(codes{i}.code{j}'), 2), 'valid');
            rf_data_decoded(:, :, 1) = temp_decoded(1:no_rf_samples, :);

            % Decode with respect to second code in current pair
            temp_decoded = conv2(rf_data_m(:, :, 2), rot90(conj(codes{i}.ccode{j}'), 2), 'valid');
            rf_data_decoded(:, :, 2) = temp_decoded(1:no_rf_samples, :);

            % Beamform image for the current line
            bf_temp = bft_beamform(Tmin, sum(rf_data_decoded, 3));
            bf_image(:, codes{i}.lineNo) = bf_image(:, codes{i}.lineNo) + bf_temp;
        end
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
imagesc([x_lines(1) x_lines(end)]*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));

title('Repeat 20 times. Correlation with Codes.');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
%colormap(gray);
caxis([-55 0]);
