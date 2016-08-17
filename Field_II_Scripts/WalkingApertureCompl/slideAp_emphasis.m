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

%% Define codes used for excitation
codes = cell(0);

% Use all same code pairs or not
% (using all the same code pair results in lots of interference)
useAllSamePairs = 0;

% Load codes (have nice cross correlation properties)
tempLoad = load('C:\Users\User\Dropbox\Grad_School\Summer Codes\OvernightPairGeneration\CompPairs_From_Nonlinear_Optimizer\compPairs_len_10_simMain.mat');
allCodes = tempLoad.('pairsSoFar');

% Used if we choose to make all the code pairs the same
firstCode =    allCodes(1,:); %[ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
secCode =   allCodes(2,:); %[ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)]; %

% Evenly space transmit focus locations in x
% numTransFocLocs = 5; % Set odd for symmetry in transmit focusing
% transFocSpacing = 6e-3; % (m)
% transFocStart = -(numTransFocLocs-1)/2*transFocSpacing;

% Set the transmit focus locations
transFocLocs = [0/1000 0 40/1000    
                 6/1000 0 40/1000
                 -6/1000 0 40/1000
                ];
numTransFocLocs = size(transFocLocs,1);
            
for i = 1:numTransFocLocs
    if(useAllSamePairs == 0)
        codes{i}.code = allCodes(i*2-1,:);
        codes{i}.ccode = allCodes(i*2,:);    
    elseif(useAllSamePairs == 1)
        codes{i}.code = firstCode;
        codes{i}.ccode = secCode;
    end
    
    %codes{currCodeInd}.focus = [transFocStart+transFocSpacing*(i-1) 0 50/1000]; % A transmit focus location (x,y,z)
    codes{i}.focus = transFocLocs(i,:);
end

%% Initial setup (ex. transducer, phantom, image range properties) 

% Store all the transmit focus locations 
transFocLocs = zeros(length(codes),3);
for i = 1:length(codes)
    transFocLocs(i,:) = codes{i}.focus;
end

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
impulse_response = sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
% pht_pos = [10/1000 0 50/1000
%            -10/1000 0 50/1000]; % Position
pht_pos = transFocLocs;
pht_amp = 20*ones(size(pht_pos,1),1); %[20; 20]; % Strength of scattering

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

% Set receive focus points
% Each row is a point in the form (x,y,z)
zDepth = 50e-3; % (m)
xFocusVals = -15e-3:1e-4:15e-3; % (m)
reFocus = zeros(numel(xFocusVals), 3);
reFocus(:,3) = zDepth;
reFocus(:,1) = xFocusVals;


% Set number of beamforming lines
no_lines = size(reFocus,1); 

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
leftIndApForTF = zeros(no_CodePairs,startLastAp);

for currTransmit = 1:2
    % Get RF data from each aperture (as the aperture slides along) 
    % Each column has data from one transducer element
    RFTransRunTot = zeros(no_rf_samples_c,no_elements);    
   
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
            xdc_apodization(rcv, 0, ones(1,no_elements));

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
            RFTransRunTot = RFTransRunTot + codeRF;

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

% Beamform at each point of interest
bfImag = zeros(no_rf_samples,no_lines); 
for currLine = 1:no_lines      
    % Set current receive focus point
    currRe = reFocus(currLine,:);   
    focusX = currRe(1);

    % Set beamforming apodization for this aperture
    %[apo_vector_tx,apo_vector_rx] = getApodization(focusX,no_active_tx,no_active_rx,width, kerf, no_elements);
    %bft_apodization(xdc, 0, apo_vector_rx);

    % Set current receive focus        
    bft_center_focus([focusX 0 0]); % Set up start of line for dynamic focus       
    bft_dynamic_focus(xdc, 0, 0);    % Set direction of dynamic focus  
    
    %for currCodeIndex = 1:no_CodePairs
        % Beamform and store resulting line
        % (beamform data decoded with respect to relevant code)        
        currLineBf = bft_beamform(Tmin, sum(sumDec, 3));     
        bfImag(:,currLine) = bfImag(:,currLine) + currLineBf;    
    %end
    
end

% Normalize and do envelope detection
env_bf = abs(hilbert(bfImag));
env_bf = env_bf / max(max(env_bf));
toPlot = 20*log10(env_bf+eps);
% Plot image
figure;
imagesc([min(reFocus(:,1)) max(reFocus(:,1))]*1000, [Rmin Rmax]*1000, toPlot);
title('Beamformed Image: Different Pairs');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);

field_end
bft_end

