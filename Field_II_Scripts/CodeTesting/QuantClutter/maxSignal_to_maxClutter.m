clear

% Load figure
% currDir = pwd; % Grab current directory
% cd('C:\Users\User\Dropbox\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\CodeTesting\AddBeamsTesting')
% [FileName,PathName] = uigetfile('*.fig','Select the MATLAB figure file');
% figLoaded = open(FileName);
%cd(currDir); % Change back to current directory

% Always load same figure for debug
figLoaded = open(['C:\Users\User\Dropbox\Grad_School\',...
'HighSpeedCodes\CompNOCodes\Field_II_Scripts\',...
'CodeTesting\AddBeamsTesting','\Tran32_ReceiveBeamAll_16ParallelX_5ParallelZ.fig']);

% Get data from loaded figure
h = gcf;
hChildren = h.Children;  
dataObjs = get(hChildren, 'Children'); %handles to low-level graphics objects in axes
imgDataObj = dataObjs(2); % Not very robust! (Would like to call by name)
imgData = imgDataObj{:}.get('CData');

%% Set region of signal and of clutter
% Center region rectangle
widthCenter = 1; % mm
heightCenter = 0.5; % mm

sigCenterMm = [0,40 + heightCenter / 2]; % (x,y) in mm

% % Side region rectangle
% widthSide= 1; % mm
% heightSide = 0.5; % mm
% 
% sideCenterMm = [0,40 + heightSide / 2]; % (x,y) in mm

%% Find corresponding region in loaded data
% Convert from mm to rows and columns
yAxis = imgDataObj{:}.get('YData');
yRange = [yAxis(1) yAxis(end)];

xAxis= imgDataObj{:}.get('XData');
xRange = [xAxis(1) xAxis(end)];

mmPerLineX = (xRange(2) - xRange(1))/size(imgData,2);
mmPerLineY = (yRange(2) - yRange(1))/size(imgData,1);
mmPerLine = [mmPerLineX mmPerLineY];

% Find position of central region in data array
mmFromLeft = (sigCenterMm(1) - 0) + (0 - xAxis(1));
elsFromLeft = mmFromLeft / mmPerLineX;

mmFromTop = sigCenterMm(2) - yAxis(1);
elsFromTop = mmFromTop / mmPerLineY;

sigCenterArray = round([elsFromLeft, elsFromTop]);

% Get data in central region
widthLines = widthCenter / mmPerLineX;
heightLines = heightCenter / mmPerLineY;%
centrData = imgData(sigCenterArray(2) +(-floor(heightLines/2):floor(heightLines/2) ),sigCenterArray(1) +(-floor(widthLines/2):floor(widthLines/2)));

% Find position of side region in data array
mmFromLeft = (sigCenterMm(1) - 0) + (0 - xAxis(1));
elsFromLeft = mmFromLeft / mmPerLineX;

mmFromTop = sigCenterMm(2) - yAxis(1);
elsFromTop = mmFromTop / mmPerLineY;

sigCenterArray = round([elsFromLeft, elsFromTop]);

% Find max value in central region
maxCent = max(max(centrData));
disp('Max in center:')
disp(maxCent)

% Debug - shows region we grabbed data in
imgData(sigCenterArray(2) +(-floor(heightLines/2):floor(heightLines/2) ),sigCenterArray(1) +(-floor(widthLines/2):floor(widthLines/2))) = 0;

figure
imagesc(xRange, yRange, imgData);
title('Selected Region');
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')
colorbar
colormap(gray);
caxis([-55 0]);


