path('/home/tjh/git/zemp_lab/Matlab/bft',path);
path('/home/tjh/Matlab Simulator',path);
load('/media/share/Tyler/Data/sa_neck.mat');

f0 = Trans.frequency*1e6;              %  Central frequency                        [Hz]
fs = f0*Receive(1).samplesPerWave;     %  Sampling frequency                       [Hz]
c = Resource.Parameters.speedOfSound;              %  Speed of sound                           [m/s]
no_elements = Trans.numelements;      %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = Trans.spacing*lambda;    % Pitch - center-to-center        [m]
width = Trans.elementWidth*lambda;     % Width of the element      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]

%  Initialize the program

bft_init;

%  Set some paramters

bft_param('c', c);
bft_param('fs', fs);


% Create some apertures.
xdc = bft_linear_array(no_elements, width, kerf);

%  Define and create the image
no_lines = 256;

imgWidth = PData(1).Size(2)*lambda;
xScanlines = linspace(-imgWidth/2, imgWidth/2, no_lines);

nSamples = 2*(Receive(1).endDepth-Receive(1).startDepth)*Receive(1).samplesPerWave;

if (exist('usTimeStart','var'))
    start_t = usTimeStart;
else
    start_t = 0;
end

bft_no_lines(no_lines);

%It turns out that this order is super important... interesting...
for j = 1 : no_lines
  bft_apodization(xdc, 0 , ones(1,no_elements),j)
  bft_center_focus([xScanlines(j) 0 0],j);
  bft_dynamic_focus(xdc, 0, 0,j);
end

%
% Make one low-resolution image at a time and sum them
%
element_x = [-(no_elements-1)/2:(no_elements-1)/2]* pitch;

tic
for emission_no = 120 : 120%no_elements
  %disp(['emission no: ' num2str(emission_no)]);  
  
  scat = double(RcvData{1}(nSamples*(emission_no-1)+1:nSamples*emission_no,:,end));

  %beamformed = bft_beamform_dsta(start_t,scat, emission_no);
  beamformed = bft_beamform_dsta(start_t,scat, [element_x(emission_no), 0, 0 ]);
  
  if (emission_no == 1)
      bf_image = beamformed;
  else
      bf_image = bf_image + beamformed;
  end
end

fprintf('bft: %f\n', toc);

tic
%{
for emission_no = 1 : no_elements
  disp(['emission no: ' num2str(emission_no)]);  

  [x,y,beamformed] = slidingDS(RcvData{1}(nSamples*(emission_no-1)+1:nSamples*emission_no,:,end));
  
  if (emission_no == 1)
      image = beamformed;
  else
      image = image + beamformed;
  end
end

figure;
imagesc(xScanlines*1000,1:size(image,1)/fs*c/2*1000,20*log10(abs(hilbert(image/max(max(image)))) + 0.001))
axis('image'); colormap(gray)
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
title('Beamformed by in-house BFT');
%}
fNum = 0;
[x,y,image] = slidingDSSA(RcvData{1}(:,:,1));
fprintf('in-house: %f\n', toc);


%  Release the allocated memory
bft_end

%
%  Display the image
%

bf_image = abs(hilbert(bf_image));                  % Envelope detection
bf_image = bf_image / max(max(bf_image));


image = abs(hilbert(image));                  % Envelope detection
image = image / max(max(image));

figure;
imagesc(xScanlines*1000,1:1000/fs*c/2*1000,20*log10(bf_image(1:1000,:) + 0.001))
axis('image'); colormap(gray)
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
title('Beamformed by BFT ');

figure;
imagesc(xScanlines*1000,1:1000/fs*c/2*1000,20*log10(image(1:1000,:) + 0.001))
axis('image'); colormap(gray)
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
title('Beamformed by in-house ');
