% Set apodization (walking aperture)(from Antares_psf_dyn_rec)
% xFocus = x coordinate of focus point
%   This script focuses at [xFocus, 0, 0]
% N_active_tx = number of transmitting transducers
% N_active_rx = number of receiving transducers

% Chooses to use the elements centered around xFocus
% Uses leftmost or rightmost elements if the focus point is too far to the
% side to use a centered apodization about xFocus

% apo_vector_tx stores how much each element is used for transmission
% apo_vector_rx stores how much each element is used for receiving
function [apo_vector_tx,apo_vector_rx] = getApodization(xFocus,N_active_tx,N_active_rx,width, kerf, N_elements)    

    % How much we activate the active transducers
    % (could add Hamming/Hanning apodization to this, if desired)
     apo_tx = ones(1,N_active_tx);
     apo_rx = ones(1,N_active_rx);

    % Number of elements before those used for transmit
    N_pre_tx = round(xFocus/(width+kerf) + N_elements/2 - N_active_tx/2);
    if(N_pre_tx < 0)
        N_pre_tx = 0;
    end;
    % Number of elemnets after those used for transmit
    N_post_tx = N_elements - N_pre_tx - N_active_tx;
    if(N_post_tx < 0)
        N_post_tx = 0;
        N_pre_tx =  N_elements - (N_post_tx + N_active_tx);
    end    
    % State how each element will be used for transmit
    apo_vector_tx = [zeros(1, N_pre_tx) apo_tx zeros(1, N_post_tx)];
    
    % Number of elements before those used for receive
    N_pre_rx = round(xFocus/(width+kerf) + N_elements/2 - N_active_rx/2);
    if(N_pre_rx < 0)
        N_pre_rx = 0;
    end;
    % Number of elements after those used for receive
    N_post_rx = N_elements - N_pre_rx - N_active_rx;
    if(N_post_rx < 0)
        N_post_rx = 0;
        N_pre_rx =  N_elements - (N_post_rx + N_active_rx);
    end 
    % State how each element will be used for receive
    apo_vector_rx = [zeros(1, N_pre_rx) apo_rx zeros(1, N_post_rx)];
end