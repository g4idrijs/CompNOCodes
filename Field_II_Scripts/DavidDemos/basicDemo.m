%% Image a simple phantom
% We want to use our codes as the impulse response
% and we want to avoid using the beamforming

% Edited by David Egolf.
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
N_elements      = 192;            % Number of physical elements
N_active_tx     = 64;             % Number of active tx elements
rx_fnum         = 2.1;            % Receive f-number

%% Set number of active receive elements for constant receive F-number
% F-number = imaging depth / aperture size
% Aperture size = N_active*sizeEachTransducer
imgDepth = rx_focus(3);
rx_ap           = imgDepth/rx_fnum;
N_active_rx     = round(rx_ap/(width+kerf)); 

%% Start the Field II program, 
field_init(0);

%% Set sampling frequency and use triangles for the aperture
set_sampling(fs);
set_field('use_triangles', 0);

%% Set emission aperture 
emit_aperture = xdc_focused_array(N_elements, width, element_height, ...
kerf, Rfocus, no_sub_x, no_sub_y, tx_focus);

%% Set impulse response and excitation of the emit aperture
% NEED TO PUT OUR CODES IN HERE

% Emitted ultrasound field at some POINT(?) if we input an impulse
% function into the transducer.
% Shouldn't this change as a function of position?

% Set impulse response to an impulse function
% impulse_response = sin(2*pi*f0*(0: 1/fs :number_cycles/f0));
lenSendExcite = 30; %length(2*pi*f0*(0: 1/fs :number_cycles/f0));
impulse_response = zeros(1,lenSendExcite);
impulse_response(1) = 1;
%figure; plot(0:length(impulse_response)-1,impulse_response, 'o');
xdc_impulse(emit_aperture, impulse_response);

% What we input to the aperture
% https://www.mathworks.com/matlabcentral/answers/46898-repeat-element-of-a-vector-n-times-without-loop
codeA1 = [0.0522411047763500,-0.0463065222031407,-0.159079283718754,-0.0693107676478534,0.0326691538381036,-0.366967215994592,-0.0386693020035286,-0.0143369692611974,0.140207954406141,-0.00676366072602007]';
codeA2 = [0.0109838703772579,-0.230561048945342,0.0882281042602866,0.0366598057588537,-0.000730468678358921,0.0788491864122134,-0.150228859791945,-0.110127580451643,-0.0196758712765472,0.0321803341860585]';



curr_line = 1; % Line to get data on
codeA1Line = walkingApImg(codeA1,lenSendExcite,emit_aperture,N_elements,width,...
    element_height,kerf,Rfocus,no_sub_x,no_sub_y,rx_focus,impulse_response,...
    N_active_tx, N_active_rx, tx_focus, curr_line);
subplot(2,2,1)
plot(codeA1Line)
title('Code A1')

codeA2Line = walkingApImg(codeA2,lenSendExcite,emit_aperture,N_elements,width,...
    element_height,kerf,Rfocus,no_sub_x,no_sub_y,rx_focus,impulse_response,...
    N_active_tx, N_active_rx, tx_focus, curr_line);
subplot(2,2,2)
plot(codeA2Line)
title('Code A2')

subplot(2,2,3)
plot(codeA1Line)
hold on
plot(codeA2Line)
legend('Code 1','Code 2')
title('Code A1 and A2 Overlaid')

subplot(2,2,4)
plot(codeA1Line + codeA2Line)
title('Sum of responses')





