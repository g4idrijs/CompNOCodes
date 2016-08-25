%PHASED_DYN_IMAGE Create a phased-array B-mode image, using the 
%    commands for setting a dynamic focusing

%VERSION 1.0, 29 Feb 2000, Svetoslav Nikolov

path('/home/tjh/git/zemp_lab/Matlab/bft',path);
path('/home/tjh/git/fieldII',path);

f0 = 5e6;              %  Central frequency                        [Hz]
fs = 20e6;            %  Sampling frequency                       [Hz]
c = 1540;              %  Speed of sound                           [m/s]
B = .6;               %  Relative bandwith                        [fraction]
no_elements = 128;      %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = lambda / 2;    % Pitch - center-to-center                  [m]
width = .95*pitch;     % Width of the element                      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]
height = 10/1000;      % Size in the Y direction                   [m]
 
 
%  Define the impulse response of the transducer
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';
excitation = impulse_response;

%  Define the phantom

pht_pos = [0 0 20;
           0 0 30;
           0 0 40;
           0 0 50;
           0 0 60;
           0 0 70;
           0 0 80;] / 1000;         %  The position of the phantom
pht_amp = 20*ones(7,1);      %  The amplitude of the back-scatter

%  Define the focus 
focus_r = [20;30;40;50;60;70;80;90] / 1000;
T = (focus_r-5/1000)/c *2;

%  Initialize the program
field_init(0);
bft_init;

%  Set some paramters
set_field('c', c);
bft_param('c', c);

set_field('fs', fs);
bft_param('fs', fs);


% Create some apertures.

xmt = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);
rcv = xdc_linear_array(no_elements,width,height,kerf,1,1,[0 0 0]);

xdc = bft_linear_array(no_elements, width, kerf);


% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

xdc_excitation(xmt, excitation);


% Set the apodization

xdc_apodization(xmt, 0, ones(1,no_elements))
xdc_apodization(rcv, 0, ones(1,no_elements))
bft_apodization(xdc, 0 , ones(1,no_elements))

%  Define and create the image
sector = 30 * pi / 180;
no_lines = 256;
d_theta = sector / (no_lines-1);
theta = -(no_lines-1) / 2 * d_theta;

Rmax = max(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) + 15/1000;

no_rf_samples = ceil(2*Rmax/c * fs);
rf_line = zeros(no_rf_samples, 1);
bf_line = zeros(no_rf_samples, 1);

env_line = zeros(no_rf_samples, no_lines);
env_bf =  zeros(no_rf_samples, no_lines);


xmt_r = (max(focus_r) + min(focus_r) )/2;
bf = cell(no_lines,1);
for i = 1 : no_lines
  rf_line(:) = 0;  
  theta;
  
  xmt_f = [sin(theta)*xmt_r, zeros(length(xmt_r),1), cos(theta)*xmt_r];
  xdc_center_focus(xmt,[0 0 0])
  xdc_center_focus(rcv,[0 0 0])
  bft_center_focus([0 0 0]);
  
  xdc_focus(xmt, 0, xmt_f);
  xdc_dynamic_focus(rcv, 0, theta, 0);
  
  %  Beamform with Field II
  [rf_temp, t(i)] = calc_scat(xmt,rcv, pht_pos, pht_amp);
  
  %  Beamform with BFT
  bft_dynamic_focus(xdc, theta, 0)
  xdc_focus_times(rcv, 0, zeros(1,no_elements));
  [rf_data, start_t] = calc_scat_multi(xmt,rcv, pht_pos, pht_amp);
  
  rf_data = [zeros(300,no_elements); rf_data; zeros(300,no_elements)];

  start_t = start_t - 300  / fs;
  bf_temp = bft_beamform(start_t, rf_data);
  
  start_sample = t(i)*fs; no_temp_samples = length(rf_temp);
  
  rf_line(start_sample:start_sample+no_temp_samples-1) = rf_temp(1:no_temp_samples);
  env_line(:,i) = abs(hilbert(rf_line(:)));

  start_sample = floor(start_t*fs); no_temp_samples = length(bf_temp);
  bf{i} = bf_temp;
  
  bf_line(start_sample:start_sample+no_temp_samples-1) = bf_temp(1:no_temp_samples);
  env_bf(:,i) = abs(hilbert(bf_line(:)));
  theta = theta + d_theta;

end

%  Release the allocated memory

field_end
bft_end
env_line = env_line / max(max(abs(env_line)));
env_bf = env_bf / max(max(abs(env_bf)));

figure;
subplot(1,2,1)
imagesc([-sector/2 sector/2]*180/pi,[0 Rmax]*1000,20*log10(env_line + 0.001))
axis('image')
xlabel('Angle [deg]');
ylabel('Axial distance [mm]')
title('Beamformed by Field II ');

subplot(1,2,2)
imagesc([-sector/2 sector/2]*180/pi,[0 Rmax]*1000,20*log10(env_bf + 0.001));
title('Beamformed by BFT');
xlabel('Angle [deg]');
ylabel('Axial distance [mm]')
axis('image')

colorbar
colormap(gray)

%clc
disp([' ' 10 10 10 10 ]);
disp([9 '*****************************************************']);
disp([9 '*                                                   *']);
disp([9 '* The image beamformed by Field II is in "env_line" *']);
disp([9 '* The image beamformed by BFT is in "env_bf"        *']);
disp([9 '*                                                   *'])
disp([9 '*****************************************************']);
disp([' ' 10 10 ]);

