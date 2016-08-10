% Antares_psf_dyn_rec.m

% This program simulates point spread functions due to pulse-echo B-mode 
% imaging using a Siemens Antares ultrasound system.

% By Roger Zemp. Updated 9 Mar 2011.

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

%pitch = lambda / 2;    % Pitch - center-to-center
%width = .95*pitch;     % Width of the element
%kerf = pitch - width;  % Inter-element spacing
%element_height = 10/1000;      % Size in the Y direction

% Receive aperture
rx_ap           = rx_focus(3)/rx_fnum;
% Number of active receive elements for constant F#
N_active_rx     = round(rx_ap/(width+kerf)); 

field_init(0);

% set the sampling frequency and use triangles
set_sampling(fs);
set_field('use_triangles', 0);

% generate aperture for emission
emit_aperture = xdc_focused_array(N_elements, width, element_height, ...
kerf, Rfocus, no_sub_x, no_sub_y, tx_focus);

% Set the impulse response and excitation of the emit aperture
impulse_response = sin(2*pi*f0*(0:1/fs:number_cycles/f0));
%figure; plot(impulse_response);
xdc_impulse(emit_aperture, impulse_response);

excitation = sin(2*pi*f0*(0:1/fs:number_cycles/f0));
xdc_excitation(emit_aperture, excitation);

% Generate aperture for reception
receive_aperture = xdc_focused_array(N_elements, width, element_height, ...
kerf, Rfocus, no_sub_x, no_sub_y, rx_focus);

% Set the impulse response for the receive aperture
xdc_impulse(receive_aperture, impulse_response);

% Load the computer phantom
[phantom_positions, phantom_amplitudes] = rz_point_phantom(5/1000, 30/1000, 5);
%[phantom_positions, phantom_amplitudes] = cyst_pht(500);

% Do linear array imaging
no_lines = 251; % Number of lines in image (odd)
image_width = 10/1000;  % size of image sector
d_x = image_width/(no_lines-1); % Increment for image

% Make the different simulations
apo_tx = ones(1, N_active_tx);
apo_rx = ones(1, N_active_rx);

z_focus = tx_focus(3); % transmit focus

x = -image_width/2;
image_data = 0;
for i = 1:no_lines
    i;
    % Set the focus for this A scan
    xdc_center_focus(emit_aperture, [x 0 0]);
    xdc_focus(emit_aperture, 0, [x 0 z_focus]);
    xdc_center_focus(receive_aperture, [x 0 0]);
    xdc_focus(receive_aperture, 0, [x 0 rx_focus(3)]);
    
    % Calculate the apodization
    N_pre_tx = round(x/(width+kerf) + N_elements/2 - N_active_tx/2);
    N_post_tx = N_elements - N_pre_tx - N_active_tx;
    apo_vector_tx = [zeros(1, N_pre_tx) apo_tx zeros(1, N_post_tx)];
    
    N_pre_rx = round(x/(width+kerf) + N_elements/2 - N_active_rx/2);
    N_post_rx = N_elements - N_pre_rx - N_active_rx;
    apo_vector_rx = [zeros(1, N_pre_rx) apo_rx zeros(1, N_post_rx)];
       
    xdc_apodization(emit_aperture, 0, apo_vector_tx);
    xdc_apodization(receive_aperture, 0, apo_vector_rx);
    
    % comment this out to turn off DRF
    xdc_dynamic_focus(receive_aperture, 0, 0, 0);
    
    % Calculate the received response
    [v t1] = calc_scat(emit_aperture, receive_aperture, phantom_positions, ...
    phantom_amplitudes);
    
    % Store the result
    image_data(1:max(size(v)), i) = v;
    times(i) = t1;
    
    % Steer in another angle/scan-line
    x = x + d_x;
end

mk_psf;
figure;
imagesc(lateral, disttxrx, log(eps+abs(hilbert(psf)))); %colormap(gray);