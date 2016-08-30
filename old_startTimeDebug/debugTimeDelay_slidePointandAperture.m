% See if we can figure out how calc_scat_multi determines start_timeCode
% (The time associated with the data coming back)
% Put a single transducer below a phantom point, and slide this all the way
% along. See if the returned time delay stays constant.

% We process only one transmit focus at a time
% (to avoid any cross correlation issues)

% Image a simple phantom 
% Use the following method: (with complementary codes: one pair per transmit focus point)
%   -Slide along a specified width aperture
%   -If a transmit focus point is inside its x range, transmit focus it
%       at that point and then fire it.
%   -Receive on the elements symmetrically around the scatterer
%   -Time adjust and add together the results from this transmit
%   -Decode the data from this transmit
%   -Beamform the data from this transmit
%   -Do the same for the second transmit, and then add the two results

%% Look at time delays as scattering point changes
clear

%  Initialize Field II and BFT
field_init(-1);
bft_init;

pointSpotsX = -8e-3:1e-3:8e-3; % Phantom/aligned aperture positions

% Hold vector of start times from calc_scat_multi
% (holds one time for each phantom / aligned aperture pair)
start_timeCodeVect = [];
figure
% Go through each phantom / aperture position
for pointCount = 1:numel(pointSpotsX)
clearvars -global -except pointSpotsX pointCount start_timeCodeVect
addpath('C:\Users\User\Dropbox\Grad_School\Summer Codes\GitCodes\Field_II_Scripts\WalkingApertureCompl')

% Current phantom position in x
phtX = pointSpotsX(pointCount);

%% Define codes and transmit focus points
codes = cell(0);

% Load codes (have nice cross correlation properties)
%tempLoad = load('C:\Users\User\Dropbox\Grad_School\Summer Codes\OvernightPairGeneration\CompPairs_From_Nonlinear_Optimizer\compPairs_len_10_simMain.mat');
%allCodes = tempLoad.('pairsSoFar');
allCodes = [repelem([1 1 1 -1],1,1);
            repelem([1 1 -1 1],1,1)];

% Set the transmit focus locations (each row is an (x,y,z) point)
zTransDepth = 40e-3; % (m)

xTransFocusVals = phtX; %-8e-3:1e-3:8e-3; %-6e-3:1e-4:6e-3; % (m)
transFocLocs = zeros(numel(xTransFocusVals), 3);
transFocLocs(:,3) = zTransDepth;
transFocLocs(:,1) = xTransFocusVals;
numTransFocLocs = size(transFocLocs,1);
  
% Currently using same code for all transmit focus points
for i = 1:numTransFocLocs
    codes{i}.code = allCodes(1,:); %allCodes(i*2-1,:);
    codes{i}.ccode = allCodes(2,:); %allCodes(i*2,:);    
    %codes{i}.focus = transFocLocs(i,:);
end

% Current phantom position (changes on each loop)
pht_pos = [phtX 0 40/1000];
                
pht_amp = 20*ones(size(pht_pos,1),1); % Strength of scattering

%% Initial setup (ex. transducer, phantom, image range properties) 

% Transducer properties
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 200;     %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]

%  Define the impulse response of the transducer
impulse_response = [1]; %sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Set Field and BFT parameters
set_field('c', c);  bft_param('c', c);
set_field('fs', fs); bft_param('fs', fs);

% Create transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Create beamforming array
xdc = bft_linear_array(no_elements, width, kerf);

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

%% Imaging preparation (ex. sliding aperture setup, receive focus, beamforming lines)

% Don't receive focus using Field II (will use beamforming toolbox)
xdc_focus_times(rcv, 0, zeros(1,no_elements));

% The number of code pairs is the number of transmit focus points
no_CodePairs = length(codes);
no_TransFocPoints = no_CodePairs;

% Set the size of the sliding aperture
wSlide = 5e-3; % width (m)
%no_elSlide = round((wSlide + kerf) / (width + kerf)); % Number of elements in sliding aperture
no_elSlide = 1;
no_active_tx = no_elSlide;
no_active_rx = no_elSlide;

%% Imaging
% The two transmits happen one after the other
%  Index of starting element of rightmost aperture
startLastAp = no_elements-no_elSlide;

% Focus on each transmit focus location
% (we process these one at a time to get a zero cross correlation image)
 no_lines = size(transFocLocs,1);

decodedRunTot = 0;
for transFocInd = 1:no_lines
    % Get current transmit focus
    transFoc = transFocLocs(transFocInd,:);
    xFocus = transFoc(1);

    % Calculate aperture directly below this transmit focus
    [apo_vector_tx,apo_vector_rx] = getApodization(xFocus,no_active_tx,no_active_rx,width, kerf, no_elements);
    
    % Set transmit and receive apodization as this aperture
    xdc_apodization(xmt, 0, apo_vector_tx);
    xdc_apodization(rcv, 0, apo_vector_rx);          

    plot(apo_vector_tx)
    title('Apodization Position')
    pause(0.05)
    
    % Transmit focus 
    xdc_focus(xmt, 0, transFoc)

    % Set excitation code (currently always using same code)
    currCodes = codes{1}; % DEBUG    
    
    % Fire each code in the complementary pair   
    for currTransmit = 1:2                   
        % Get RF scattering data by firing...                
        % ...first code in pair
        if(currTransmit == 1)            
            %currEl = find(apo_vector_tx);
            %waveform = [1];
            %ele_waveform (xmt, currEl, waveform);
            
            xdc_excitation(xmt, currCodes.code);
            [codeRF, start_timeCode] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);             
            
            start_timeCodeVect = [start_timeCodeVect start_timeCode];
        else
            %...second code in pair
            % xdc_excitation(xmt, currCodes.ccode);
            [codeRF, start_timeCode] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp); 
        end     
    end    
end 

end

% Plot time delays as a function of position along aperture
figure
plot(pointSpotsX*1000,start_timeCodeVect )
xlabel('Center of aperture and phantom X position (mm)')
ylabel('Calc\_scat\_multi start time (s)')
title('Start time for Aligned Aperture As Aperture and Phantom Point Slide')

% Release memory
field_end
bft_end

