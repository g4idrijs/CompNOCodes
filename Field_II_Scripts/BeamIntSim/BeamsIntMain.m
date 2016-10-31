% Specify two subapertures
% Fire on one and receive on the other
% (No codes yet)

% I wonder if it's possibly not focusing on the subaperture properly
% but instead using the entire array for apodization...
% Can avoid this issue by using only one element per subaperture

% Need to use exponential decay

clear
close all
addpath('../Field_II')

%% Subaperture properties
numeElSubAp = 100;
numSubApCentPad = 10; % How many apertures to put between the two interfering ones

%% Attenuation properties
% (using settins from user guide, p. 16)
set_field ('att',1.5*100);
set_field ('Freq_att',0.5*100/1e6);
set_field ('att_f0',3e6);
set_field ('use_att',1);

%% Transducer properties
fs = 1*100e6;           % Simulation sampling frequency [Hz]
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
c = 1540;               % Speed of sound                           [m/s]
no_elements = numeElSubAp*2 + numeElSubAp*numSubApCentPad;    % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 

no_active_tx = numeElSubAp;     % Number of active elements for transmit 
                                % sub-aperture       
                                
Ts = 1 /fs; % Sampling period
%Note: depending on Ts and number_cycles/f0, impulse_response may not have complete cycles
impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0)); 
impulse_response = impulse_response.*hanning(length(impulse_response))';

%% Set phantom
pht_pos = [0 0 28; 0 0 30; 0 0 32; 0 0 34];
pht_amp = ones(size(pht_pos,1),1);

%% Set up Field II
field_init(-1);
set_field('c', c);
set_field('fs', fs);

%% Set Field II transmit and receive arrays
% We will use apodization to send on one and receive on the other
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Set impulse response
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

% Set apodization
for mismatch = [0, 1] % Receive and send on different subapertures
    if(mismatch == 1)
        xdc_apodization(xmt, 0, [ones(1,numeElSubAp) zeros(1,numeElSubAp*numSubApCentPad) zeros(1,numeElSubAp)]) % Transmit on left subaperture
        xdc_apodization(rcv, 0, [zeros(1,numeElSubAp) zeros(1,numeElSubAp*numSubApCentPad) ones(1,numeElSubAp)]) % Receive on right subaperture
    else
        xdc_apodization(xmt, 0, [ones(1,numeElSubAp) zeros(1,numeElSubAp*numSubApCentPad) zeros(1,numeElSubAp)]) % Transmit on left subaperture
        xdc_apodization(rcv, 0, [ones(1,numeElSubAp) zeros(1,numeElSubAp*numSubApCentPad) zeros(1,numeElSubAp) ]) % Receive on right subaperture
    end

    % Set center point for focusing [NOT SURE ABOUT THIS]
    % (for calculating delay times, and as starting point for dynamic focusing)
    % Each subarray focuses in a different spot
    padWidth = numSubApCentPad*width;
    xdc_center_focus(xmt,[-padWidth/2 - width/2 0 0])
    xdc_center_focus(rcv,[padWidth/2 + width/2 0 0])

    % Turn on dynamic receive focusing
    % xdc_dynamic_focus(rcv, 0, 0, 0);

    %% Carry out simulation
    [scat, start_time] = calc_scat(xmt, rcv, pht_pos, pht_amp);

%     scat = abs(hilbert(scat));
    if(mismatch == 0)
        maxNoMismatch = max(max(abs(scat)));
    end
    scat = scat / maxNoMismatch;
    if(mismatch == 0)
        scatMatch = scat;
    else
        scatMismatch = scat;
    end
    plot(scat);
    hold on
    
    title('Signal Received from Adjacent Beam')
    ylabel('Normalized Intensity')
    xlabel('Time Index')
end

legend('No mismatch','Mismatched')
disp('Ratio of maxima (mismatch / no mismatch):')
ratioMax = max(max(scatMismatch))/ max(max(scatMatch)) - 1










