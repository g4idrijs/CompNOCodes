% Use complementary codes to image a simple phantom.
% Weird sidelobes
% -Likely caused by using calc_scat
% --We probably need to be doing our own beamforming
% --Why would beamforming before cross correlating mess up things?

% Modified by David Egolf.
% close all

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
fs              = 200e6;          % Samping frequency [Hz]
c               = 1540;           % Speed of sounds [m/s]
lambda          = c/f0;           % Wavelength [m]
width           = 0.2/1000;       % Width of element [m]
element_height  = 0.05/1000;         % Height of element [m]
kerf            = 0.02/1000;      % Kerf [m]
tx_focus        = [0 0 30]/1000;  % Transmit focus [m]
Rfocus          = 30/1000;        % Elevation focal radius
rx_focus        = [0 0 30]/1000;  % Receive focus [m]
no_sub_x        = 1;              % Number of x subdivisions of element
no_sub_y        = 15;             % Number of y subdivisions of element
N_elements      = 192;              % Number of transducer array elements
N_active_tx     = 64;             % Number of active tx elements
rx_fnum         = 2.1;            % Receive f-number

%% Start the Field II program, 
field_init(0);

%% Set sampling frequency
set_sampling(fs);
% Using triangles means exact solutions are used (I believe - FAQ)
set_field('use_triangles', 0);

%% Set emit aperture 
% emit_aperture = xdc_linear_array(N_elements, width, element_height, ...
% kerf, no_sub_x, no_sub_y, tx_focus);

% Try a focused array (to eliminate grating lobes as a concern)
emit_aperture = xdc_focused_array(N_elements, width, element_height, ...
  kerf, Rfocus, no_sub_x, no_sub_y, tx_focus);

%% Set emit aperature impulse response
% Set impulse response... 

% ...to be a sine wave
% tVect = linspace(0,2*pi,30);
% impulse_response = sin(tVect);

% ...to be an impulse function (so codes are sent un-modulated: simple)
impulse_response = [ones(1,32)]; 

% Plot impulse response (what does length correspond to?)
% figure; plot(0:length(impulse_response)-1,impulse_response, 'o');

% Feed the produced response to the aperature
xdc_impulse(emit_aperture, impulse_response);

%% Set receive aperture
% receive_aperture = xdc_linear_array(N_elements, width, element_height, ...
% kerf, no_sub_x, no_sub_y, rx_focus);
receive_aperture = xdc_focused_array(N_elements, width, element_height, ...
  kerf, Rfocus, no_sub_x, no_sub_y, tx_focus);

%% Set receive aperature impulse response (using transmit impulse response)
xdc_impulse(receive_aperture, impulse_response);

%% Load the computer phantom
% calling form: rz_point_phantom(dz, z_start, Npoints)
[phantom_positions, phantom_amplitudes] = rz_point_phantom(5/1000, 30/1000, 5);

%% Loop through focus spots (get one line for each)
xFocusSpots = [-0.014:1e-4:0.014];
lineMat = cell(numel(xFocusSpots),1);

% Store apodization values
% Number of receive transducers (for a desired F number)
rx_ap           = rx_focus(3)/rx_fnum;
N_active_rx     = round(rx_ap/(width+kerf)); 
if(N_active_rx > N_elements)
    N_active_rx = N_elements;
end

for j = 1:numel(xFocusSpots)

    %% Set focus (at a single point)
    xFocus = xFocusSpots(j); % Where in the x plane do we take this image?

    % Set up emit aperature focus
    xdc_center_focus(emit_aperture, [xFocus 0 0]);     
    xdc_focus(emit_aperture, 0, [xFocus 0 tx_focus(3)]);
    % Set up transmit aperature focus
    xdc_center_focus(receive_aperture, [xFocus 0 0]); 
    xdc_focus(receive_aperture, 0, [xFocus 0 rx_focus(3)])

    % We send and receive on all of the elements
%     apo_vector_tx = ones(1,N_elements);
%     apo_vector_rx = ones(1,N_elements);
%     xdc_apodization(emit_aperture, 0, apo_vector_tx);
%     xdc_apodization(receive_aperture, 0, apo_vector_rx);

    % Set apodization (walking aperture)(from Antares_psf_dyn_rec)

    % How much we activate the active transducers
     apo_tx = ones(1,N_active_tx);
     apo_rx = ones(1,N_active_rx);

    % The number of elements before those that are activated
    % Seems to have issues if we want to focus to the side...
    % (N_pre_tx < 0)
    % Uses the elements centered around xFocus
    % Uses leftmost or rightmost elements if  the focus point is to the
    % side of the transducer    
    N_pre_tx = round(xFocus/(width+kerf) + N_elements/2 - N_active_tx/2);
    if(N_pre_tx < 0)
        N_pre_tx = 0;
    end;
    N_post_tx = N_elements - N_pre_tx - N_active_tx;
    if(N_post_tx < 0)
        N_post_tx = 0;
        N_pre_tx =  N_elements - (N_post_tx + N_active_tx);
    end    
    apo_vector_tx = [zeros(1, N_pre_tx) apo_tx zeros(1, N_post_tx)];
    
    N_pre_rx = round(xFocus/(width+kerf) + N_elements/2 - N_active_rx/2);
    if(N_pre_rx < 0)
        N_pre_rx = 0;
    end;
    N_post_rx = N_elements - N_pre_rx - N_active_rx;
    if(N_post_rx < 0)
        N_post_rx = 0;
        N_pre_rx =  N_elements - (N_post_rx + N_active_rx);
    end 
    apo_vector_rx = [zeros(1, N_pre_rx) apo_rx zeros(1, N_post_rx)];
       
    xdc_apodization(emit_aperture, 0, apo_vector_tx);
    xdc_apodization(receive_aperture, 0, apo_vector_rx);


%     % Try a Hamming apodization (fire on all, but weighted)
%     apo_vector_tx = hamming(N_elements)';
%     apo_vector_rx = hamming(N_elements)';
%     xdc_apodization(emit_aperture, 0, apo_vector_tx);
%     xdc_apodization(receive_aperture, 0, apo_vector_rx);   
    

    %% Specify codes
    codeA1 = [ones(1,8), ones(1,8), ones(1,8), -ones(1,8)]; % 1 1 1 -1 
    codeA2 = [ones(1,8), ones(1,8), -ones(1,8), ones(1,8)]; %  1 1 -1 1
%     codeA1 = [1]; codeA2 = [1];
    codes = [codeA1; codeA2];

    %% Image with each code   
    ccfMat = []; % Store cross correlation data
    
    for i = 1:size(codes,1)           
        %% Set current emit aperture excitation    
        currCode = codes(i,:);
        xdc_excitation(emit_aperture, currCode);      
        
        % Set dynamic focusing on receive
        xdc_dynamic_focus(receive_aperture, 0, 0, 0);
        
        % Calculate received voltage
        % This beamforms the data (time adjusts, weights, and sums)
        [rf_data, t1] = calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
                
        % Calculate cross correlation with currrent code
        ccfMat(i,:) = xcorr(rf_data,currCode);    

    end
    %% Sum results of cross correlation with each code  
    lineMat{j} = sum(ccfMat); % Store the current line
end

% Pad shorter sequence with zeroes at the end
maxLen = 0;

% Find the length of the longest sequence
for i = 1:size(lineMat,1)
    currLen = size(lineMat{i},2);
    if(currLen > maxLen || maxLen == 0)
        maxLen = currLen;
    end
end

% Pad to max length
for i = 1:size(lineMat,1)
   currLen = size(lineMat{i},2); 
   toPad = maxLen - currLen;
   lineMat{i} = [lineMat{i}, zeros(1,toPad)];
end

% Arrange the data in columns
imgData = cell2mat(lineMat)';

% Plot the data
figure
toPlot = imgData/max(max(imgData));
toPlot = 20*log10(abs(toPlot) + eps);
imagesc(toPlot)
colormap(gray);
colorbar
caxis([-40 0]);