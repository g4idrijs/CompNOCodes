% Finds signal received on one subaperture when a different one fires
% Both subapertures are focused axially around their centers
% Based on SingleElEnerg.m

clear

addpath('../Field_II')

%% Transducer properties
fs = 1*100e6;           % Simulation sampling frequency [Hz]
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
c = 1540;               % Speed of sound                           [m/s]
no_elements = 192;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 
                                 
Ts = 1 /fs; % Sampling period

impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0)); 
impulse_response = impulse_response.*hanning(length(impulse_response))';

%% Set up Field II, including attenuation
field_init(-1);
set_field('c', c);
set_field('fs', fs);

% From Field II user guide (some plausible default values?)
% Set the attenuation to 1.5 dB/cm and 0.5 dB/[MHz cm] around 3 MHz
set_field ('att',1.5*100);
set_field ('Freq_att',0.5*100/1e6);
set_field ('att_f0',3e6);
set_field ('use_att',1);


%% Create arrays
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Set impulse response
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

%% Set up subaperture properties

% Size of subapertures (number of elements)
numElSubAp = 64;

% Choose transmit subaperture leftmost element index
% (left = small, right = large)
tranSubInd = 1;

% Get indices of transmit elements
transElemsInd = tranSubInd:(tranSubInd + numElSubAp - 1);

% Calculate x positions of transmit subaperture
xmt_info = xdc_get(xmt, 'rect');
xElements = xmt_info(8, :);
xPosTransAp = xElements(transElemsInd);

%% Set up transmit focus
% Find x coordinate of center of transmitting subaperture
xCentTrans = xPosTransAp(round(numElSubAp/2));

% Set transmit focus point
zFocusTrans = 30e-3; % m
transFocus = [xCentTrans 0 zFocusTrans];

% Set phantom at transmit focus point
pht_pos = transFocus;
pht_amp = ones(size(pht_pos,1),1);

%% Set up transmit focus excitation
% Set excitation pattern
excit = impulse_response;

% Will store excitation for each element
% (one columns for each element's excitation over time, which is row-wise)
forceMaxDelay = 5+389+2+3; 
max_code_length = numel(excit);
beam = zeros(forceMaxDelay+max_code_length, no_elements);

% Get actual voltages to feed to transmit subarray
beam(:, transElemsInd) = beam(:, transElemsInd) + focusBeam(excit, transFocus, xPosTransAp, fs, c, forceMaxDelay);

%% Carry out transmit focus excitation
% Set up excitation
ele_waveform(xmt, (1:no_elements)', beam(:, :)');

% Fire the beam and record the scattering
[rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);


%% Choose receive apodization
% Keep track of maximum data at different apodizations
maxDataVect = [];

rightMostInd = no_elements - numElSubAp;    % Left-most index of rightmost subarray
apodIndCollection = 1:rightMostInd;         % Apodization left-most indices to loop through
for rcvSubInd = apodIndCollection           % Loop through different apodizations
    rcvElemsInd =  rcvSubInd:(rcvSubInd + numElSubAp - 1);  % Choose elements to receive on
    apoData = rf_data(:,rcvElemsInd);                       % Get data that was received by those elements
    maxData = max(max(apoData));                            % Find the maximum voltage recorded in this apodization
    
    maxDataVect = [maxDataVect maxData];
end

% Plot received data, if only one apodization explored
if(numel(apodIndCollection) == 1)
    figure
    plot(apoData);
    title(sprintf('Transmit on Subaperture %d, receive on Subaperture %d',tranSubInd, rcvSubInd))
    xlabel('Time index')
    ylabel('Voltage Magnitude')
else % Plot maximum voltages otherwise
    figure
    plot(apodIndCollection, maxDataVect/max(max(abs(maxDataVect))))
    ylabel('Maximum Energy Received (normalized) [no receive focusing]')
    xlabel('Index of Leftmost Element of Subaperture')
    title('Energy Leakage due to Beam Transmitted by Leftmost Transducer')
end

%% Get magnitude of maximum voltage
maxVolt = max(max(rf_data));

%%
% %% Align RF data 
% % Minimum and maximum distance to image
% Rmax = 35e-3;
% Rmin = 15e-3;
% 
% % Minimum and maximum time to wait while imaging
% Tmin = 2*Rmin / c; Tmax = 2*Rmax / c;
% 
% % Minimum and maximum samples to take
% % (if we started collecting at time t = 0)
% Smin = floor(Tmin * fs); Smax = ceil(Tmax * fs);
% 
% % Align RF data in data
% no_rf_samples = Smax - Smin + 1;
% alignedRF = alignRF(rf_data,start_time,fs,Smin,Smax,no_rf_samples,no_elements);
% %% Beamform RF data







