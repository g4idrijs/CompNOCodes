% Transducer properties
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
fs = 100e6;             % Sampling frequency                       [Hz]
c = 1540;               % Speed of sound                           [m/s]
no_elements = 192;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 

%  Define the impulse response of the transducer
Ts = 1 /fs; % Sampling period
impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

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

xdc_excitation(xmt, impulse_response);

xdc_focus(xmt, 0, [0 0 40/1000]);

% Get excitation at a range of positions over time
xRange = [-5:0.1:5]/1000; % m
zRange = [35:0.1:45]/1000; % m    
pointsInterest = zeros(numel(zRange),3);
pointsInterest(:,3) = zRange'; 

% Determine indices of samples to keep
[Rmax, Rmin, ~, ~, ~, Smin_cVis, Smax_cVis, ~, no_rf_samples_cVis] =...
calcSampleTimeRanges(pointsInterest, codes, c, fs);
Smin_cVis = round(Smin_cVis/2); % Adjust for no echo here
Smax_cVis = round(Smax_cVis/2);
no_rf_samples_cVis = round(no_rf_samples_cVis /2);

for xRangeInd = 1:numel(xRange)
   for zRangeInd = 1:numel(zRange)
       currX = xRange(xRangeInd);
       currZ = zRange(zRangeInd);
       [currResp,start_timeVis] = calc_hp(xmt,[currX 0 currZ]);              
       plot(currResp) 
       currResp = alignRF(currResp,start_timeVis,fs,Smin_cVis,Smax_cVis,no_rf_samples_cVis,1);

       % Store sequence at right (x,z) location
       respMat(zRangeInd,xRangeInd,:) = currResp;
   end
end

plot(currResp)   

%% Show the animation
numFrames = size(respMat,3);
for frame = 300:numFrames 
   imagesc([min(xRange), max(xRange)]*1000,[Rmin, Rmax]*1000, respMat(:,:,frame));

   title(['First Code Set. Line: ', num2str(lineNo),'. Frame: ', num2str(frame), ' of ', num2str(numFrames)]);
   colorbar;
   caxis([-1 1]*1e-11)
   ylabel('z (mm)')
   xlabel('x (mm)')

   pause(0.01);
end