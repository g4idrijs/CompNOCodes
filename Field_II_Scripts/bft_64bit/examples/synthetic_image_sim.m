path('/home/tjh/git/zemp_lab/Matlab/bft',path);
path('/home/tjh/git/fieldII',path);

f0 = 4e6;              %  Central frequency
fs = 100e6;            %  Sampling frequency
c = 1540;              %  Speed of sound
no_elements = 65;      %  Number of elements in the transducer

lambda = c / f0;       % Wavelength
pitch = lambda / 2;    % Pitch - center-to-center
width = .95*pitch;     % Width of the element
kerf = pitch - width;  % Inter-element spacing
height = 10/1000;      % Size in the Y direction
 
 
%  Define the impulse response of the transducer
impulse_response = sin(3*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(length(impulse_response))';
excitation = sin(2*pi*f0*(0:1/fs:3/f0));
excitation = excitation .* hanning(length(excitation))';
%  Define the phantom

%pht_pos = [12 0 40; 0 0 40] / 1000;   %  The position of the phantom
pht_pos = [0 0 40] / 1000;   %  The position of the phantom

[m n] = size(pht_pos);
pht_amp = 20*ones(m,1);      %  The amplitude of the back-scatter

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

no_lines=256;
d_x_line = width*no_elements / (no_lines-1);
x_line = -(no_lines-1) / 2 * d_x_line;

% Set the impulse responses
xdc_impulse(rcv, impulse_response);
xdc_impulse(xmt, impulse_response);

xdc_excitation(xmt, excitation);


%  Set the delays for one whole image
%
bft_no_lines(no_lines);

for i = 1 : no_lines
  %bft_apodization(xdc,0,hanning(no_elements)',i);
  bft_apodization(xdc, 0 , ones(1,no_elements),i);
  bft_center_focus([x_line 0 0],i);
  bft_dynamic_focus(xdc, 0, 0,i);
  x_line = x_line + d_x_line;
end
%
%  Allocate memory for the image
%
%Rmax = max(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) + 40/1000;
%Rmin = min(sqrt(pht_pos(:,1).^2 + pht_pos(:,2).^2  + pht_pos(:,3).^2)) - 10/1000;
%if (Rmin < 0) Rmin = 0; end;
%Tmin = 2*Rmin / c; Tmax = 2*Rmax / c;
%Smin = floor(Tmin * fs); Smax = ceil(Tmax * fs);

%no_rf_samples = Smax - Smin + 1;

%bf_image = zeros(no_rf_samples, no_lines);
%
% Make one low-resolution image at a time and sum them
%

xdc_focus_times(xmt,0,zeros(1,no_elements)); 
xdc_focus_times(rcv,0,zeros(1,no_elements));


element_x = [-(no_elements-1)/2:(no_elements-1)/2]* pitch;

Smin=1;
Smax = 6000;

tic
for emission_no = 1 : 65
  disp(['emission no: ' num2str(emission_no)]);  
  xdc_apodization(xmt,0,[zeros(1,emission_no-1) 1 zeros(1, no_elements - emission_no)]);

  [scat, start_time] = calc_scat_multi (xmt, rcv, pht_pos, pht_amp);

  start_sample = floor(start_time * fs + 0.5);
  end_sample = start_sample + max(size(scat))-1;
  
  scat = [zeros(start_sample - Smin, no_elements); scat; zeros(Smax - end_sample,no_elements)];
 
  start_time = 0;
  %beamformed = bft_beamform_dsta(start_time,scat, emission_no);
  beamformed = bft_beamform_dsta(start_time,scat, [element_x(emission_no), 0, 0 ]);
  
  if (emission_no==1)
      bf_image = beamformed;
  else
      bf_image = bf_image + beamformed;
  end
end

t = toc



%  Release the allocated memory

field_end
bft_end



%
%  Dispplay the image
%

bf_image = abs(hilbert(bf_image));                  % Envelope detection
bf_image = bf_image / max(max(bf_image));

figure;
imagesc(1000*(1:no_lines)*d_x_line,1000*(start_time+(1:size(bf_image,1))/fs*c)/2,20*log10(bf_image + 0.001))
axis('image')
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
title('Beamformed by BFT ');

clc
disp([' ' 10 10 10 10 ]);
disp([9 '*****************************************************']);
disp([9 '*                                                   *']);
disp([9 '* The image beamformed by BFT is in "bf_image"      *']);
disp([9 '*                                                   *'])
disp([9 '*****************************************************']);
disp([' ' 10 10 ]);
