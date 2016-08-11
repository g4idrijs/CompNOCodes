% FLASH_COMPLEMENTARY_SIM 
% Simulate a flash image with complementary codes at many focal zones
% VERSION 1.0, August 10, 2016 (Tarek)
%
% This version uses BFT to beamform images
% 
% Version 1.1, August 11, 2016 (David)
% Adds calcuation of axial and lateral full width half maximum
% as well as signal to noise ratio

%% Simulation properties (initial setup) 
clear

% Add white Gaussian noise
useNoise = 1;

% Define codes used for excitation
codes = cell(0);
codes{1}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
codes{1}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];

% tempDat = load('bestPairsSQP41.mat');
% x = tempDat.('x');
% codes{1}.code = x(1,:); 
% codes{1}.ccode = x(2,:);

codes{1}.focus = [0 0 50/1000]; 

no_lines = 257;

% Transducer properties
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]

rx_fnum         = 0; % 2.1;   % Receive f-number - Set to 0 to turn off
                       % constants F-number reconstruction

%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:1/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';

%  Define the phantom
pht_pos = [0 0 50/1000]; % Position
pht_amp = 20; % Strength of scattering

% Calculate minimum and maximum samples and times for phantom
Rmax = max(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) + 5/1000;
Rmin = min(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) - 5/1000;
if (Rmin < 0) Rmin = 0; end;

Tmin = 2*Rmin / c; Tmax = 2*Rmax / c;
Smin = floor(Tmin * fs); Smax = ceil(Tmax * fs);
max_code_length = max(cellfun(@(x) max(length(x.code), length(x.ccode)), codes));
Smin_c = Smin; Smax_c = Smax + max_code_length + 1000;

no_rf_samples = Smax - Smin + 1;
no_rf_samples_c = Smax_c - Smin_c + 1;

%  Initialize Field II and BFT
field_init(-1);
bft_init;

%  Set Field and BFT parameters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);

bft_no_lines(no_lines);

% Create transmitting and receiving apertures
xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

% Create beamforming array
xdc = bft_linear_array(no_elements, width, kerf);

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

%% Set the apodization
xdc_apodization(xmt, 0, ones(1,no_elements))
xdc_apodization(rcv, 0, ones(1,no_elements))

%% Simulate transmit and receive

% Set a center reference point for focusing
xdc_center_focus(xmt,[0 0 0])
xdc_center_focus(rcv,[0 0 0])

% Don't receive focus
xdc_focus_times(rcv, 0, zeros(1,no_elements));

% Set up the properties of each line we beamform along
d_x_line = width*no_elements / (no_lines-1); % Line width
x_line = -(no_lines-1) / 2 * d_x_line; % Position of leftmost line
z_points = linspace(Rmin, Rmax, 100); % For constant F number
no_points = length(z_points);
T = z_points/c*2;

% Set up beamforming
for i = 1 : no_lines
   % Set reference point for dynamic focusing for this line
   % Dynamic focusing focuses everywhere in this line
    bft_center_focus([x_line 0 0], i); 
    
    % Set angle of dynamic focus line for this line
    % We set it to point straight at the focus
    bft_dynamic_focus(xdc, 0, 0, i);  
    
    % Constant F-Num reconstruction
    if (rx_fnum == 0)
        % Use all elements if fnum = 0
        bft_apodization(xdc, 0, ones(1, no_elements), i);
    else
        rx_ap           = z_points/rx_fnum;
        no_active_rx    = min(round(rx_ap/(width+kerf)), no_elements);

        mid_element = round((x_line+d_x_line*(no_lines-1)/2)/(width)+0.5);
        start_element = round(max(mid_element - no_active_rx/2 + 1, 1));
        end_element = round(min(mid_element+no_active_rx/2, no_elements));

        apo_vector_rx = zeros(no_points, no_elements);
        for j = 1:no_points
            apo_vector_rx(j, start_element(j):end_element(j)) = 1;
        end
    
        bft_apodization(xdc, T, apo_vector_rx, i);
    end
  
    x_line = x_line + d_x_line;
end

%% Simulate imaging with each code
% rf_data generated from imaging with each code and its associated focus
% zone is added together to simulate imaging with all codes. This is then 
% repeated for complementary versions.

% Rf data matrix holds non-complementary and complementary simulated
% RF data.
rf_data_m = zeros(no_rf_samples_c, no_elements, 2);
for i = 1:2*length(codes) % Go through both first code and second code for each pair
        
    % Extract code and focus information for current code
    if (i > length(codes)) % Working with second code "ccode"
        code = codes{i-length(codes)}.ccode;
        focus = codes{i-length(codes)}.focus;
        col = 2; % Where to store results
    else
        code = codes{i}.code; % Working with first code "code"
        focus = codes{i}.focus;
        col = 1; % Where to store results
    end        
    
    xdc_focus(xmt, 0, focus);
    xdc_excitation(xmt, code);
    [rf_data, start_time] = calc_scat_multi(xmt, rcv, pht_pos, pht_amp); 
    
    start_sample = floor(start_time * fs)+1;
    end_sample = start_sample + size(rf_data, 1) - 1;

    if (start_sample < Smin_c)
      rf_data = rf_data(Smin_c-start_sample+1:end, :);
      start_sample = Smin_c;
    end
    rf_data = [zeros(start_sample - Smin_c, no_elements); rf_data; zeros(Smax_c - end_sample, no_elements)];

    if size(rf_data, 1) > no_rf_samples_c
      rf_data = rf_data(1:no_rf_samples_c, :);
    end

    rf_data_m(1:size(rf_data, 1), :, col) = rf_data_m(1:size(rf_data, 1), :, col) + rf_data;
end

% Add noise
if(useNoise == 1)
    % Add white gaussian noise (prior to beamforming)     
    noisePowerdBW =  -4.45e2; % Use -inf to turn off noise
    rf_data_m(:,:,1)  = rf_data_m(:,:,1) + wgn(size(rf_data_m,1), size(rf_data_m,2),noisePowerdBW);
    rf_data_m(:,:,2)  = rf_data_m(:,:,2) + wgn(size(rf_data_m,1), size(rf_data_m,2),noisePowerdBW);
end

%% Decode images
% For each code, cross-correlate rf_data_sum with each code and add the
% result to rf_data_decoded
rf_data_nc_decoded = 0;
rf_data_c_decoded = 0;
for i = 1:length(codes)
    %temp_decoded = xcorr2(rf_data_m(:, :, 1), codes{i}.code');
    %temp_decoded = temp_decoded(length(codes{i}.code):end, :);
    temp_decoded = conv2(rf_data_m(:, :, 1), rot90(conj(codes{i}.code'), 2), 'valid');
    temp_decoded = temp_decoded(1:no_rf_samples, :);
    rf_data_nc_decoded = rf_data_nc_decoded + temp_decoded;

    %temp_decoded = xcorr2(rf_data_m(:, :, 2), codes{i}.ccode');
    %temp_decoded = temp_decoded(length(codes{i}.ccode):end, :);
    temp_decoded = conv2(rf_data_m(:, :, 2), rot90(conj(codes{i}.ccode'), 2), 'valid');
    temp_decoded = temp_decoded(1:no_rf_samples, :);
    rf_data_c_decoded = rf_data_c_decoded + temp_decoded;
end

% Can either add decoded data for non-complementary and complementary then
% beamform.
rf_data_decoded = rf_data_nc_decoded + rf_data_c_decoded;
bf_temp = bft_beamform(Tmin, rf_data_decoded);

% Or we can beamform decoded data then add results
%bf_temp = bft_beamform(Tmin, rf_data_nc_decoded);
%bf_temp = bf_temp + bft_beamform(Tmin, rf_data_c_decoded);

%% Calculate signal to noise ratio
if(useNoise == 1)
    % Signal / sd(noise)
    % std(noise) =  about 6e-21 with Gaussian white noise at -4.45e2 dBW.
    % and with the following codes:
    %     codes{1}.code = [ones(1, 8), ones(1, 8), ones(1, 8), -ones(1, 8)];
    %     codes{1}.ccode = [ones(1, 8), ones(1, 8), -ones(1, 8), ones(1, 8)];
    % (see the findNoise script)
    maxSig = max(max(bf_temp));
    stdNoise = 6e-21;
    sigToNoise = maxSig/stdNoise;
    disp('Signal to noise ratio:');
    disp(sigToNoise);
end

% Perform envelope detection and normalize
env_bf = abs(hilbert(bf_temp));
env_bf = env_bf / max(max(env_bf));

%% Calculate full width half maximum resolution
% (Needs generalization to work with a several point phantom)
% Find the index of the maximum signal
[maxVal,maxLocLin] = max(env_bf(:));
[maxLocZ, maxLocX] = ind2sub(size(env_bf),maxLocLin);

% Axial (z direction)
zIntVect = env_bf(:,maxLocX);
fwhmIndicesZ = fullWidthHM(zIntVect); % FWHM in indices
lenPerIndexZ = (Rmax - Rmin)/numel(zIntVect);
fwhmLengthZ = lenPerIndexZ*fwhmIndicesZ; % Convert to length
disp('Axial full width half maximum: (mm)')
disp(fwhmLengthZ*1000)

% Lateral (x direction)
xIntVect = env_bf(maxLocZ,:);
fwhmIndicesX = fullWidthHM(xIntVect); % FWHM in indices
lenPerIndexX = (no_lines*d_x_line)/numel(xIntVect);
fwhmLengthX = lenPerIndexX*fwhmIndicesX; % Convert to length
disp('Lateral full width half maximum: (mm)')
disp(fwhmLengthX*1000)

%% Plot image
figure;
imagesc([-1/2 1/2]*no_lines*d_x_line*1000, [Rmin Rmax]*1000, 20*log10(env_bf+eps));
title('Beamformed Image');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray);
caxis([-55 0]);
