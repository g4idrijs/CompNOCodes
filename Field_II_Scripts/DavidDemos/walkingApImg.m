% Returns the data in the no_lines line, currently, summed but not
% beamformed (I think)

% NEEDS EDITING

function sumDataCurr = walkingApImg(transCode,lenSendExcite,emit_aperture,N_elements,width,...
    element_height,kerf,Rfocus,no_sub_x,no_sub_y,rx_focus,impulse_response,...
    N_active_tx, N_active_rx, tx_focus,no_lines)


excitation = repmat(transCode,1,lenSendExcite/numel(transCode))';
excitation = excitation(:)';
%excitation = sin(2*pi*f0*(0:1/fs:number_cycles/f0));

% seqSend = [1,-1];
% baseShape = ones(1,lenSendExcite/numel(seqSend));
% excitation = [];
% for i = 1:length(seqSend)
%     excitation = [excitation seqSend(i)*baseShape];
% end

%excitation = [ones(1,lenSendExcite/2) -ones(1,lenSendExcite/2)];
xdc_excitation(emit_aperture, excitation);

%% Generate aperture for reception
% Use focused array 
% (curved in z direction so as to focus)
receive_aperture = xdc_focused_array(N_elements, width, element_height, ...
kerf, Rfocus, no_sub_x, no_sub_y, rx_focus);

%% Set the impulse response for the receive aperture
% Currently is same as for transmit aperture
xdc_impulse(receive_aperture, impulse_response);

%% Load the computer phantom
% rz_point_phantom(dz, z_start, Npoints)
[phantom_positions, phantom_amplitudes] = rz_point_phantom(5/1000, 30/1000, 5);

%% Set up linear array imaging
% x direction = direction linear array extends in
image_width = 10/1000;  % Range of transmit focus in x direction
d_x = image_width/(no_lines-1); % Increment between lines

% Setup apodization for transmit and receive
apo_tx = ones(1, N_active_tx);
apo_rx = ones(1, N_active_rx);

% Set transmit focus location in z direction
% z direction is in the direction of depth into the phantom
z_focus = tx_focus(3); 

% x = 0 at the center of the array, and we start imaging at the left
% side of our image
x = -image_width/2;
image_data = 0;
for i = 1:no_lines % Get image data along each line    
    %% Set the focus for this A scan (this line)
    xdc_center_focus(emit_aperture, [x 0 0]); % Emit aperature
    % Sets start of focus line
    
    xdc_focus(emit_aperture, 0, [x 0 z_focus]);
    % Creates focus line
    
    xdc_center_focus(receive_aperture, [x 0 0]); % Receive aperature
    xdc_focus(receive_aperture, 0, [x 0 rx_focus(3)]);
    
    %% Calculate apodization (which elements transmit, receive)
    N_pre_tx = round(x/(width+kerf) + N_elements/2 - N_active_tx/2);
    N_post_tx = N_elements - N_pre_tx - N_active_tx;
    apo_vector_tx = [zeros(1, N_pre_tx) apo_tx zeros(1, N_post_tx)];
    
    N_pre_rx = round(x/(width+kerf) + N_elements/2 - N_active_rx/2);
    N_post_rx = N_elements - N_pre_rx - N_active_rx;
    apo_vector_rx = [zeros(1, N_pre_rx) apo_rx zeros(1, N_post_rx)];
    
    xdc_apodization(emit_aperture, 0, apo_vector_tx);
    xdc_apodization(receive_aperture, 0, apo_vector_rx);
    
    %% DRF = dynamic receive focus
    xdc_dynamic_focus(receive_aperture, 0, 0, 0);
    
    %% Calculate the received response
    
    % Carries out beamforming (take into account different delays
    % from different parts of aperture)
    % and sums image data from each part of aperture
    % Use calc_scat_multi to avoid beamforming / summing
    
% [sumDataCurr t1] = calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
%     plot(v)
%     title('Calc\_scat')    
    
    % Each column is the reponse of an element in the receive aperature
    [v t1] = calc_scat_multi(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    
    % Debug - show the raw data collected by some element
%     currEl = 63;
%     currData = v(:,currEl);
%     figure
%     plot(currData);
%     
    % Debug - sum the responses (add columns together)
    sumDataCurr = sum(v,2);
%     figure
%     plot(sumDataCurr);  
%     title('Calc\_scat\_multi')    
      
    % Store the result
    % image_data(1:max(size(v)), i) = v;
    % times(i) = t1;
    
    % Change the center of focus, as well as current apodization
    x = x + d_x;
end

%% Plot the results
% mk_psf;
% figure;
% % Plots a matrix
% % Taking log(hilbert) returns an envelope of the response
% imagesc(lateral, disttxrx, log(eps+abs(hilbert(psf)))); colormap(gray);