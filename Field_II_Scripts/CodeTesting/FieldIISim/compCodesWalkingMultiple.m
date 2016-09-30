% compCodesWalkingMultiple.m

% Authors: Tarek, David

% Description Uses Field II and BFT to simulate a walking aperture of 
% parallel focal zones in both X and Z directions using complementary 
% codes. Number of parallel zones in X direction determine how many lines 
% will be simultaneously walked, and number of parallel zones in Z 
% direction determine how many parallel focal zones in the Z direction.
clear

addpath('../../Field_II', '../../bft_64bit', '../QuantClutter');

%% Simulation properties - Higher Level

onlySimCentLine = 0; % Only simulate center line
doneCent = 0; % Have we simulated the center line already?

% Plot center line
plotCent = 1;

% Plot expected center line
plotExpCent = 0;

% Create a video of the transmitted beam (first code) in the center line
visBeam = 0; % Visualize the beam (either single spot, or animation)
plotSingLoc = 1; % Only visualize the energy at one point in the beam

 % Visualize echo energy
visEcho = 0;
zRange = 0.1/1000; %(35:-0.1:0.1)/1000;
visLocEcho = zeros(numel(zRange),3);
visLocEcho(:,3) = zRange';

% Visualize spatial impulse response
visImpRe = 0;
impPos = [0 0 40]/1000;

% Use single pair of codes for all transmissions
useCustomCode = 0;

custCode1 =  [repelem(1,50) repelem(-1,50)]; % sin(linspace(0,2*pi,100)); % triangle(2,30); % % [repelem(1,10) repelem(-1,10)]; % % First code in pair
custCode1 = custCode1.*hanning(length(custCode1))';
custCode2 = custCode1;  % Second code in pair

% Envelope the received response
useEnv = 1;

% Transducer sampling frequency [Hz]
fs = 1*100e6;             

% Number of focal zones
numCodesX = 3; % Number of parallel focal zones in X
numCodesZ = 3; % Number of parallel focal zones in Z
numFocZones = numCodesX * numCodesZ;

% Use two codes that get along well as neighbors to construct code set
useNeigh = 1; 

% Sets focal zone spacing
image_width = 12/1000;  % Image width to simulate [m]
no_lines = 151;         % Number of lines in image 

% Repeat elements in codes @numRepeat times
numRepeat = 1;

%% Simulation properties - Lower Level

% Load codes used for excitation

% Use two codes that are good neighbors
if(useNeigh == 1)
%     tempLoad =  load('../../../Complementary Pairs/NeighborCodes_Sept29_69neigh.mat');
    tempLoad =  load('../../../Complementary Pairs/len10_2codes_1150.mat');
    codeSet = tempLoad.('x');
    if(size(codeSet,1) ~= 2*2)
       error('Only two neighbors allowed for now.') 
    end
    
    codeSet =  codeFromNeigh(codeSet(1:2,:), codeSet(3:4,:), numCodesX, numCodesZ);
else

%     tempLoad =  load('../../../Complementary Pairs/len10_16codes_minInt2.mat');
%     tempLoad =  load('../../../Complementary Pairs/len100_10codes.mat');
    tempLoad =  load('../../../Complementary Pairs/len10_9codes3P3.mat');
    codeSet = tempLoad.('x');
end

lenCodes = size(codeSet,2);

% Store codes in a cell array for later use
codes = cell(0);
codesToUse = [1:2:numCodesX*numCodesZ*2-1]; % Which codes to use from code set

if(numCodesZ == 1)
   focalPoints_z = 40/1000; 
else
   focalPoints_z = linspace(30, 50, numCodesZ)/1000;
end
 
% Store codes
for i = 1:numCodesX
    for j = 1:numCodesZ   
        if(useCustomCode == 1) % Use user specified codes as desired (for all pairs)
            codes{i}.code{j} = custCode1;
            codes{i}.ccode{j} = custCode2;
        else
            codes{i}.code{j} = repelem(codeSet(codesToUse(1+(i-1)*numCodesZ+(j-1)), :),1,numRepeat); % 1st pair code
            codes{i}.ccode{j} =  repelem(codeSet(codesToUse(1+(i-1)*numCodesZ+(j-1))+1, :),1,numRepeat); % 2nd pair code
        end
        
        codes{i}.focusZ(j) = focalPoints_z(j);
    end
end

% Transducer properties
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
c = 1540;               % Speed of sound                           [m/s]
no_elements = 192;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 

no_active_tx = 64;     % Number of active elements for transmit 
                        % sub-aperture
rx_fnum_constant = 0;  % F-number for constant f-num reconstruction.
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

forceMaxDelay = 5+389+2+3;    % Forces the maximum delay in manual focusing to 
                        % this value. This is needed to align multiple 
                        % focused beams in the Z direction. [samples]

%  Define the impulse response of the transducer
Ts = 1 /fs; % Sampling period
%Note: depending on Ts and number_cycles/f0, impulse_response may not have complete cycles
impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0)); 
impulse_response = impulse_response.*hanning(length(impulse_response))';
if(isequal(impulse_response,1))
   usedDelta = 1; 
else
    usedDelta = 0;
end
                        
%% Phantom definition

%  Define the phantompht
pHWidth = 4; % Phantom half width (mm)
pht_pos = [ 
           0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
           -pHWidth 0 30; -pHWidth 0 35; -pHWidth 0 40; -pHWidth 0 45; -pHWidth 0 50;
           pHWidth 0 30; pHWidth 0 35; pHWidth 0 40; pHWidth 0 45; pHWidth 0 50;
           ] / 1000; %  The position of the phantom points (x,y,z)
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

% Create Field II transmitting and receiving apertures
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
    % Jump straight to middle line if desired
    if (onlySimCentLine == 1)
        if(doneCent == 0)
            lineNo = round((numel(x_lines)/length(codes))/2); 
            doneCent = 1;
        else
            break
        end
    end
    
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
    
    % Don't transmit or receive focus using Field II
    xdc_focus_times(xmt, 0, zeros(1, no_elements));
    xdc_focus_times(rcv, 0, zeros(1, no_elements));
    
    % Set transducer excitation and simulate with calc_scat_multi
    % Align scatter data (force data to lie in a specified time window)
    % The corresponding sample window is [Smin_c, Smax_c]
    % Also need to shift rf_data by forceMaxDelay before or after aligning 
    % to account for focusing delays.
    ele_waveform(xmt, (1:no_elements)', beam(:, :, 1)'); 
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);
    rf_data = [rf_data(1+forceMaxDelay:end, :); zeros(forceMaxDelay, size(rf_data, 2))];
    rf_data = alignRF(rf_data,start_time,fs,Smin_c,Smax_c,no_rf_samples_c,no_elements);
    rf_data_m(:, :, 1) = rf_data;
    
    % Visualize the beam from the first code
    if(visBeam == 1 && lineNo == round(focalZoneSpacing_x/2))   
        % animInfo stores animation information
        % allows us to grab different frame later
        animInfo = visBeamSubModule(xmt,codes,c,fs,plotSingLoc);
        doneVisBeam = 1;
    end
    
    % Plot echo energy in the middle line as it passes through visLocEcho
    if(visEcho == 1 && lineNo == round(numel(x_lines)/2))
        visEchoSubModule(visLocEcho, xmt, rcv);
    end
    
    % Visualize spatial impulse response
    if(visImpRe == 1 && lineNo == round(focalZoneSpacing_x/2))
       [h, ~] = calc_h(xmt, impPos);
       figure
       plot(h)
       title(sprintf('Spatial impulse response at depth of %d mm', (impPos(3)*1000)))
       ylabel('Spatial impulse response (m/s)')
       xlabel('Time Index (right is later in time)')       
    end
    
    
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

% Perform envelope detection and normalize
if(useEnv == 1)
    env_bf = abs(hilbert(bf_image));
else
    env_bf = abs(bf_image);
end

env_bf = env_bf / max(max(env_bf));

%% Plot central line
if(plotCent == 1)
    figure
    midLine = env_bf(:,round(numel(x_lines)/2)); 
    midLine = midLine/max(max(midLine));
    xAxisSimLine = linspace(Rmin,Rmax, numel(midLine))*1000;
    plot(xAxisSimLine, 20*log10(midLine+eps));

    if(usedDelta == 1)
        title(sprintf('Central Line of PSF using Delta Transducer Impulse Response \n %d Focal Zone(s). Code length: %d. Times repeated: %d. f_s: %d MHz.', numFocZones, lenCodes, numRepeat, fs/(1e6)))
    else
        title(sprintf('Central Line of PSF with Image Width %.0f mm, %.0f Lines. \n %d Focal Zone(s). Code length: %d. Times repeated: %d. f_s: %d MHz. \n Neighbors: %i', image_width*1000, no_lines, numFocZones, lenCodes, numRepeat, fs/(1e6),useNeigh))
    end

    xlabel('Axial Distance (mm)')
    ylabel('Normalized Intensity (dB)')    

    hold on
end
%% Show expected central line
if(plotExpCent == 1)
    % Get transmit impulse response at scatterer location
    [spatialImpRespTrans,startTimeImp] = calc_h(xmt, [0 0 40]/1000);
    spatialImpRespTrans = spatialImpRespTrans(1:50);        
    spatialImpResp = [1];
    
    % Get receive impulse response at scatterer location
    [spatialImpRespRcv,~] = calc_h(rcv, [0 0 40]/1000);
    spatialImpRespRcv = spatialImpRespRcv(1:50);    
    spatialImpRespRcv = [1];
    
    % Get version of codes we correlate against by convolving with:
    % transducer impulse response
    % transmit spatial impulse response
    % receive spatial impulse response
    % transducer impulse response
    travelledCode1 = conv(impulse_response,conv(spatialImpRespRcv,conv(spatialImpRespTrans,conv(impulse_response,codes{1}.code{1}))));
    travelledCode2 = conv(impulse_response,conv(spatialImpRespRcv,conv(spatialImpRespTrans,conv(impulse_response,codes{1}.ccode{1}))));
    
    % Correlate received data with sent out code
    centSlice = abs(xcorr(codes{1}.code{1},travelledCode1) + xcorr(codes{1}.ccode{1},travelledCode2));        
    
    % Normalize
    centSlice = centSlice/max(max(centSlice));  
%   
%    % Shift to start at right time, and 
%    % pad with zeroes to make same length as received data
%    zeroBefore = round(startTimeImp*fs);
%    centSlice = [zeros(1,zeroBefore), centSlice];
%    centSlice = [centSlice zeros(1,numel(midLine) - numel(centSlice))];
%    
%    plot(linspace(Rmin,Rmax, numel(centSlice))*1000,20*log10(centSlice + eps))   
%    
%    legend('calc\_scat\_multi simulation', 'Using model')
    figure
    plot(20*log10(centSlice+eps))
    title('Expected Central Line if No Cross Beam Interference')
    xlabel('Time Index ->')
    ylabel('Magnitude')
    
%     % Show autocorrelation of codes
%     figure
%     plot( 20*log10(abs(xcorr(codes{1}.code{1}) + xcorr(codes{1}.ccode{1})))+eps);
%     title('Code Autocorrelation')
%     xlabel('Index')
%     ylabel('Value (dB)')
end

%% Plot entire image
if(onlySimCentLine == 0)
    figure;
    imagesc([x_lines(1) x_lines(end)]*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));
    title(sprintf('Simulated part of PSF. Neighbors: %i. \n %d Focal Zone(s). Code length: %d. Times repeated: %d. f_s: %d MHz.',useNeigh, numFocZones,lenCodes, numRepeat, fs/(1e6)));
    xlabel('Lateral distance [mm]');
    ylabel('Axial distance [mm]')
    axis('image')

    colorbar
    %colormap(gray);
    caxis([-55 0]);
end


% Release memory used by Field II and BFT
field_end
bft_end