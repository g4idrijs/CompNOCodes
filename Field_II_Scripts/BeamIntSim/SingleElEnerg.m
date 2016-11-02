% Created by David, Tarek November 1, 2016
% Shows energy received on one element sent by another elements

clear
close all
addpath('../Field_II')

%% Transducer properties
fs = 1*100e6;           % Simulation sampling frequency [Hz]
f0 = 6.67e6;            % Central frequency                        [Hz]
number_cycles = 2;      % Number of cycles for impulse response
c = 1540;               % Speed of sound                           [m/s]
no_elements = 128;      % Number of elements in the transducer array 

width = 0.2/1000;       % Width of element [m]
height = 5/1000;        % Height of element [m]
kerf = 0.02/1000;       % Kerf [m] 
 
                                
Ts = 1 /fs; % Sampling period

impulse_response = sin(2*pi*f0*(0:Ts:number_cycles/f0)); 
impulse_response = impulse_response.*hanning(length(impulse_response))';

%% Set phantom
pht_pos = [ 0 0 30]/1000;
pht_amp = ones(size(pht_pos,1),1);

%% Set up Field II
field_init(-1);
set_field('c', c);
set_field('fs', fs);

%% Create arrays
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Set impulse response
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

%% Set apodization

% Just transmitting on one element
txApoVect =  zeros(1,no_elements);
txApoVect(1) = 1; % Element to transmit on
xdc_apodization(xmt, 0, txApoVect)

% Just receiving on one elmenet
rcvApoVect = zeros(1,no_elements);
rcvApoVect(end) = 1; % Element to receive on
xdc_apodization(rcv, 0, rcvApoVect)

xdc_center_focus(xmt,[0 0 0])
xdc_center_focus(rcv,[0 0 0])

%% Set transmit excitation
xdc_excitation(xmt,impulse_response);

%% Set dummy focus to avoid strange default values
xdc_focus(xmt,0,[0 0 10e8]);

%% Carry out simulation
[scat, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp);

rcvData = scat(:,end);
plot(rcvData);

























