% Use complementary codes to image a simple phantom.

% Modified by David Egolf.
close all

%% Original Info
% Antares_psf_dyn_rec.m

% This program simulates point spread functions due to pulse-echo B-mode 
% imaging using a Siemens Antares ultrasound system.

% By Roger Zemp. Updated 9 Mar 2011.

%% Include Field II directories
addpath('C:\Users\User\Dropbox\Grad_School\Summer Codes\GitCodes\Field_II_Scripts')
addpath('C:\Users\User\Dropbox\Grad_School\Summer Codes\GitCodes\Field_II_Scripts\Field_II')

%% Set Field II parameters
f0              = 6.67e6;         % Transducer center frequency [Hz]
number_cycles   = 2;              % number of cycles to transmit
fs              = 100e6;          % Samping frequency [Hz]
c               = 1540;           % Speed of sounds [m/s]
lambda          = c/f0;           % Wavelength [m]
width           = 0.2/1000;       % Width of element [m]
element_height  = 5/1000;         % Height of element [m]
kerf            = 0.02/1000;      % Kerf [m]
tx_focus        = [0 0 30]/1000;  % Transmit focus [m]
Rfocus          = 30/1000;        % Elevation focal radius
rx_focus        = [0 0 30]/1000;  % Receive focus [m]
no_sub_x        = 1;              % Number of x subdivisions of element
no_sub_y        = 15;             % Number of y subdivisions of element
%N_elements      = 192;            % Number of physical elements
N_active_tx     = 64;             % Number of active tx elements
rx_fnum         = 2.1;            % Receive f-number

%% Start the Field II program, 
field_init(0);

%% Set sampling frequency
set_sampling(fs);
% Using triangles means exact solutions are used (I believe - FAQ)
set_field('use_triangles', 0);

%% Set emit aperture 
emit_aperture = xdc_linear_array(N_elements, width, element_height, ...
kerf, no_sub_x, no_sub_y, tx_focus);

%% Set emit aperature impulse response
% Set impulse response... 

% ...to be a sine wave
% tVect = linspace(0,2*pi,30);
% impulse_response = sin(tVect);

% ...to be an impulse function (so codes are sent un-modulated: simple)
impulse_response = [1, zeros(1,29)];

% Plot impulse response (what does length correspond to?)
% figure; plot(0:length(impulse_response)-1,impulse_response, 'o');

% Feed the produced response to the aperature
xdc_impulse(emit_aperture, impulse_response);

%% Set receive aperture
N_elements = 1;
receive_aperture = xdc_linear_array(N_elements, width, element_height, ...
kerf, no_sub_x, no_sub_y, rx_focus);

%% Set receive aperature impulse response (using transmit impulse response)
xdc_impulse(receive_aperture, impulse_response);

%% Load the computer phantom
% rz_point_phantom(dz, z_start, Npoints)
[phantom_positions, phantom_amplitudes] = rz_point_phantom(5/1000, 30/1000, 5);

%% Set focus (at a single point)
xFocus = 0; % Image lined up with the center of the imaging array

% Set up emit aperature focus
xdc_center_focus(emit_aperture, [xFocus 0 0]);     
xdc_focus(emit_aperture, 0, [xFocus 0 tx_focus(3)]);
% Set up transmit aperature focus
xdc_center_focus(receive_aperture, [xFocus 0 0]); 
xdc_focus(receive_aperture, 0, [xFocus 0 rx_focus(3)])

% We send and receive on the one element
apo_vector_tx = 1;
apo_vector_rx = 1;
xdc_apodization(emit_aperture, 0, apo_vector_tx);
xdc_apodization(receive_aperture, 0, apo_vector_rx);

%% Specify codes
codeA1 = [ones(1,15), -1*ones(1,15)];
codeA2 = [ones(1,15), ones(1,15)];
codes = [codeA1; codeA2];

%% Image with each code
rfMat = [];
ccfMat = [];
figure
for i = 1:size(codes,1)           
    %% Set current emit aperature excitation    
    currCode = codes(i,:);
    xdc_excitation(emit_aperture, currCode);

    % Calculate received voltage
    % This may be time adjusted - not sure
    [rf_data, t1] = calc_scat_multi(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    rfMat(i,:) = rf_data;
    
    % Plot RF data
    subplot(3,2,1+2*(i-1))
    plot(rf_data)
    title(sprintf('RF %i',i))
    
    % Calculate cross correlation with code
    % Should cross correlate with response to code
    subplot(3,2,2+2*(i-1))

    ccfMat(i,:) = xcorr(rf_data,currCode);
    plot(ccfMat(i,:));
    title(sprintf('CCF %i',i)) 
end
%% Sum results of cross correlation with each code
subplot(3,2,6)
plot(sum(ccfMat))
title('Sum of CCF')


