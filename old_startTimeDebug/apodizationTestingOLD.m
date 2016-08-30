clear

% Load cleaned up RF decoded reflection from just one point
% (as read by some subaperture)
RFforTest = load('C:\Users\User\Dropbox\Grad_School\Summer Codes\GitCodes\Zero Interference Windows\decodedRFOneAperture');
RFforTest = RFforTest.('decodedRF');
imagesc(RFforTest)

bft_init

pht_pos = [0/1000 0 40/1000
           -6/1000 0 40/1000]; % Position

%% Define codes used for excitation
codes = cell(0);

% Load codes (have nice cross correlation properties)
tempLoad = load('C:\Users\User\Dropbox\Grad_School\Summer Codes\OvernightPairGeneration\CompPairs_From_Nonlinear_Optimizer\compPairs_len_10_simMain.mat');
allCodes = tempLoad.('pairsSoFar');
    
numFocusPoints = 2;
for i = 1:numFocusPoints
    codes{i}.code = allCodes(i*2-1,:);
    codes{i}.ccode = allCodes(i*2,:);
end
       
%% Aperture setup
f0 = 5e6;              %  Central frequency                        [Hz]
fs = 100e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
no_elements = 192;     %  Number of elements in the transducer     
lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]   

bft_param('c', c);
bft_param('fs', fs);

%% Calculate minimum and maximum samples and times for phantom
[Rmax, Rmin, Tmin, Smin, max_code_length, Smin_c,Smax_c, no_rf_samples, no_rf_samples_c] =...
    calcSampleTimeRanges(pht_pos, codes,c,fs);

%% Beamform
transmitRcvAp = zeros(1,no_elements);
transmitRcvAp(26:58) = 1;

xdc = bft_linear_array(no_elements, width, kerf);


bft_apodization(xdc, 0, transmitRcvAp);
bft_center_focus([-6/1000 0 0]);
bft_dynamic_focus(xdc, 0, 0);

currLineBf = bft_beamform(Tmin, RFforTest);
plot(currLineBf)



bft_end


