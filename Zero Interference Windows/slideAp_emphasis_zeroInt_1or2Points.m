% Demo that we can avoid cross correlation by stopping our listening
% -Proof of concept for zero interference window stuff

% Image a simple phantom 
% Use the following method: (with complementary codes: one pair per transmit focus point)
%   -Slide along a specified width aperture
%   -If a transmit focus point is inside its x range, transmit focus it
%       at that point and then fire it.
%   -Time adjust and add together the results from this transmit
%   -Decode the data from this transmit
%   -Beamform the data from this transmit
%   -Do the same for the second transmit, and then add the two results

clear

addpath('C:\Users\User\Dropbox\Grad_School\Summer Codes\GitCodes\Field_II_Scripts\WalkingApertureCompl')

%% Define codes used for excitation
codes = cell(0);

% Load codes (have nice cross correlation properties)
%tempLoad = load('C:\Users\User\Dropbox\Grad_School\Summer Codes\OvernightPairGeneration\CompPairs_From_Nonlinear_Optimizer\compPairs_len_10_simMain.mat');
%allCodes = tempLoad.('pairsSoFar');
allCodes = [repelem([1 1 1 -1],1,1);
            repelem([1 1 -1 1],1,1)];

% Set the transmit focus locations
transFocLocs = [ ...%0/1000 0 40/1000
                 -3/1000 0 40/1000
                ];
numTransFocLocs = size(transFocLocs,1);
            
for i = 1:numTransFocLocs
    codes{i}.code = allCodes(i*2-1,:);
    codes{i}.ccode = allCodes(i*2,:); 
    
    %codes{currCodeInd}.focus = [transFocStart+transFocSpacing*(i-1) 0 50/1000]; % A transmit focus location (x,y,z)
    codes{i}.focus = transFocLocs(i,:);
end

%% Initial setup (ex. transducer, phantom, image range properties) 

% Transducer properties
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     
% no_active_tx = 64;     %  Number of transmitting elements in walking aperture
% no_active_rx = 64;     %  Number of receiving elements in walking aperture (used for beamforming apodization)
                       %    We actually receive on all elements

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]

%  Define the impulse response of the transducer
impulse_response = [1]; %sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
% pht_pos = [10/1000 0 50/1000
%            -10/1000 0 50/1000]; % Position
pht_pos = [ %0/1000 0 40/1000
            -3/1000 0 40/1000
                ];
pht_amp = 20*ones(size(pht_pos,1),1); %[20; 20]; % Strength of scattering

% Calculate minimum and maximum samples and times for phantom
 [Rmax, Rmin, Tmin, Smin, max_code_length, Smin_c,Smax_c, no_rf_samples, no_rf_samples_c] =...
    calcSampleTimeRanges(pht_pos, codes,c,fs);

%  Initialize Field II and BFT
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

%% Imaging preparation (ex. sliding aperture setup, receive focus, beamforming lines)

% Don't receive focus (will use beamforming toolbox)
xdc_focus_times(rcv, 0, zeros(1,no_elements));

% The number of code pairs is the number of transmit focus points
no_CodePairs = length(codes);
no_TransFocPoints = no_CodePairs;

% Set the size of the sliding aperture
wSlide = 5e-3; % width (m)
no_elSlide = round((wSlide + kerf) / (width + kerf)); % Number of elements in sliding aperture
no_active_tx = no_elSlide;
no_active_rx = no_elSlide;

%% Imaging
% The two transmits happen one after the other
% We store and process the data from each individually
firstCodeCorr = zeros(no_rf_samples, no_elements, no_CodePairs); 
secCodeCorr = zeros(no_rf_samples, no_elements, no_CodePairs);  

%  Index of starting element of rightmost aperture
startLastAp = no_elements-no_elSlide;

% Store what apertures fired at what transmit focus points
% Row = which transmit focus point
% (fill row with indices of leftmost index in apertures that fired at that
% point)

numAp = 1; % Debug: fire with only one aperture

leftIndApForTF = zeros(no_CodePairs,startLastAp);

for currTransmit = 1:2
    numApDone = 0;
    
    % Get RF data from each aperture (as the aperture slides along) 
    % Each column has data from one transducer element
    RFTransRunTot = zeros(no_rf_samples_c,no_elements);    
    
    % Only fire on some apertures (DEBUG)
    if (numApDone >= numAp)
       break 
    end
    
    for elStartFocus = 1:startLastAp % Slide the start of the aperture along
        % x position of leftmost edge of leftmost element
        leftmostX = -(no_elements*(width+kerf))/2;
        
        % Check for a focus spot in the x range of this aperture
        xStart = leftmostX + (elStartFocus-1)*(width+kerf);
        xEnd = xStart + (no_elSlide*(width+kerf) + width); % Rightmost edge, rightmost element in current aperture
        
        % Set transmit focus point
        transFocInd = find((transFocLocs(:,1) >= xStart) .* (transFocLocs(:,1) <= xEnd));
        transFoc = transFocLocs(transFocInd,:);
              
        if (numel(transFoc) > 3)
            error('Too many places to focus!')            
        elseif(numel(transFoc) == 3) % There is one place to transmit focus            

            % Set transmit apodization for current aperture            
            apo_vector_tx = [zeros(1,elStartFocus) ones(1,no_elSlide) zeros(1,no_elements-(elStartFocus+no_elSlide))];
            
            xdc_apodization(xmt, 0, apo_vector_tx);
            %xdc_apodization(rcv, 0, ones(1,no_elements));
            % Restrict aperture to same section
            %xdc_apodization(rcv, 0, apo_vector_tx);
            % Set receive aperture as that directly below point
            [~,apo_vector_rx] = getApodization(transFoc(1),no_active_tx,no_active_rx,width, kerf, no_elements)    ;
            xdc_apodization(rcv, 0, apo_vector_rx);

            % Get the complementary code pair for this transmit focus spot
            for i = 1:no_CodePairs
                if(codes{i}.focus == transFoc)
                    currCodes = codes{i};
                    
                    % Record the leftmost element used for this transmit focus point
                    % Row = which transmit focus point
                    % Column = which leftmost element was used
                    leftIndApForTF(i,elStartFocus) = 1;
                    break;
                end
            end            

            % Set transmit focus
            xMid = (xStart + xEnd)/2;
            xdc_center_focus(xmt,[xMid, 0, 0]) % Middle of current aperture
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
            % (in each transmit we fire from all sub-apertures)
            RFTransRunTot = RFTransRunTot + codeRF;
            
            numApDone = numApDone + 1;
            % Only fire on some apertures (DEBUG)
            if (numApDone >= numAp)
               break 
            end
        end
    end
   
 
    % We have collected the data for one transmit
    % Decode to get desired data for each code pair
    for currCodePair = 1:no_CodePairs
        % Get current code for decoding this code pair
        if (currTransmit == 1)
            currDecodeCode = codes{currCodePair}.code;
        else
            currDecodeCode = codes{currCodePair}.ccode;   
        end 
        % Cross correlate
        % RFTransRunTot should be trimmed to a relevant subset for decoding
        % this code
        currPairDec = conv2(RFTransRunTot, rot90(conj(currDecodeCode'), 2), 'valid');
        currPairDec = currPairDec(1:no_rf_samples, :);
        
        % Store results
        if (currTransmit == 1)
            firstCodeCorr(:,:,currCodePair) = currPairDec;
        else
            secCodeCorr(:,:,currCodePair) = currPairDec;            
        end
    end      
end
 
% Sum the results of cross correlating with the different codes
% (2D: decoded data received by each transducer, 3rdD: which code is
% decoded with respect to)
sumDec =  firstCodeCorr + secCodeCorr; % first transmit + second transmit results

decodedRF = sum(sumDec,3);
decodedRF(680:end,:) = 0;
figure
imagesc(decodedRF)

% Beamform at each point of interest

% Set receive focus points
% Each row is a point in the form (x,y,z)
zDepth = 40e-3; % (m)
xFocusVals = -12e-3:1e-4:0; % (m)
reFocus = zeros(numel(xFocusVals), 3);
reFocus(:,3) = zDepth;
reFocus(:,1) = xFocusVals;

no_lines = size(reFocus,1); 

bfImag = zeros(no_rf_samples,no_lines);

for currLine = 1:no_lines      
    % Set current receive focus point
    currRe = reFocus(currLine,:);   
    focusX = currRe(1);

    % Set beamforming apodization for this aperture
%     [apo_vector_tx,apo_vector_rx] = getApodization(focusX,no_active_tx,no_active_rx,width, kerf, no_elements);
%     bft_apodization(xdc, 0, apo_vector_rx);
    bft_apodization(xdc, 0, apo_vector_rx); % DEBUG
    
    % Set current receive focus        
    bft_center_focus([focusX 0 0]); % Set up start of line for dynamic focus       
    bft_dynamic_focus(xdc, 0, 0);    % Set direction of dynamic focus  
    
    %for currCodeIndex = 1:no_CodePairs
        % Beamform and store resulting line
        % (beamform data decoded with respect to relevant code)        
        currLineBf = bft_beamform(Tmin, decodedRF);     
        bfImag(:,currLine) = bfImag(:,currLine) + currLineBf;    
    %end
    
end

% Normalize and do envelope detection
env_bf = abs(hilbert(bfImag));
env_bf = env_bf / max(max(env_bf));
toPlot = 20*log10(env_bf+eps);

figure
% subplot(1,2,1)
% plot(linspace(Rmin,Rmax,numel(toPlot))*1000,toPlot)
% title('Intensity of Beamformed Center')
% ylabel('Intensity (dB)')
% xlabel('Depth (mm)')
% xlim([Rmin, Rmax]*1000)
% 
% % Plot image
% subplot(1,2,2)
imagesc([min(reFocus(:,1)) max(reFocus(:,1))]*1000, [Rmin Rmax]*1000, toPlot);
title('Beamformed Center');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);
field_end
bft_end

