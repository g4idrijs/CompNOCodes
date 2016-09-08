clear

% Load figure
% currDir = pwd; % Grab current directory
% cd('C:\Users\User\Dropbox\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\CodeTesting\AddBeamsTesting')
% [FileName,PathName] = uigetfile('*.fig','Select the MATLAB figure file');
% figLoaded = open(FileName);
%cd(currDir); % Change back to current directory

% Always load same figure for debug
% figLoaded = open(['C:\Users\User\Dropbox\Grad_School\',...
% 'HighSpeedCodes\CompNOCodes\Field_II_Scripts\',...
% 'CodeTesting\AddBeamsTesting','\Tran32_ReceiveBeamAll_16ParallelX_5ParallelZ.fig']);

% Get data from loaded figure
h = gcf;
hChildren = h.Children;  
dataObjs = get(hChildren, 'Children'); %handles to low-level graphics objects in axes
imgDataObj = dataObjs(2); % Not very robust! (Would like to call by name)


%% List points to investigate
analysisPoints  = [           
           0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
           -6 0 30; -6 0 35; -6 0 40; -6 0 45; -6 0 50;
           6 0 30; 6 0 35; 6 0 40; 6 0 45; 6 0 50;
           ]; % (x,y,z) (mm)
analysisPoints = [analysisPoints(:,1) analysisPoints(:,3)];

% Will store metric results
% One result for each point
metricArr = zeros(size(analysisPoints,1),1);

for pointCount = 1:size(analysisPoints,1)
    %% Set region of signal and of clutter
    center = analysisPoints(pointCount,:); % (x,z) in mm

    % Signal box
    sigWidth = 1.5; % mm
    sigHeight = 0.5; % mm
    sigCentMm = center + [0 sigHeight/2]; % (x,z) in mm

    % Clutter box
    clutterWidth = 3; % mm
    clutterHeight = 5; % mm
    clutterCentMm = center + [0 sigHeight/2]; % (x,z) in mm

    %% Find corresponding data in loaded array (figure)

    % Get data in signal center box
    % Also get the row and column indices of the box
    [arrDataSig,boxRowIndSig,boxColIndSig] = mmBoxToArray(imgDataObj, sigCentMm, sigWidth, sigHeight);
    [X,Z] = meshgrid(boxRowIndSig,boxColIndSig);
    sigPoints = [X(:) Z(:)];

    % Get data in  clutter box
    % Also get the row and column indices of the box
    [arrDataClutter,boxRowIndClutter,boxColIndClutter] = mmBoxToArray(imgDataObj, clutterCentMm, clutterWidth, clutterHeight);
    [X,Z] = meshgrid(boxRowIndClutter,boxColIndClutter);
    clutterPoints = [X(:) Z(:)];

    % Get data in clutter box but not in signal center box
    imgData = 10.^(imgDataObj{:}.get('CData')./20); % Load entire image, removes log compression
    % Grab points in clutter, not in signal
    jCluttLoc = setdiff(clutterPoints,sigPoints,'rows'); 
    jCluttLinInd = sub2ind(size(imgData),jCluttLoc(:,1),jCluttLoc(:,2)); 
    justClutterData = imgData(jCluttLinInd);    

    %% Find max value in signal and clutter regions
    % normResultsTo = 
    %[metric, message ] = maxToMax( arrDataSig, justClutterData); 
     [metric, message ] = maxToAvg( arrDataSig, justClutterData); 
   
    
    metricArr(pointCount) = metric;
end

% Visualize the metrics handed out as an array
% We want to plot this in positions corresponding to the analysis point
% positions
yPointsChunk = 5; % mm
xPointsChunk = 6; % mm

xPoints = analysisPoints(:,1);
yPoints = analysisPoints(:,2);

xChunkInd = xPoints./xPointsChunk;
xChunkInd = xChunkInd - min(xChunkInd) + 1;

yChunkInd = yPoints./yPointsChunk;
yChunkInd = yChunkInd - min(yChunkInd) + 1;

sizePlot = [numel(unique(yPoints)), numel(unique(xPoints))];
lineIndPlot = sub2ind(sizePlot,yChunkInd,xChunkInd);

plotArr = reshape(metricArr(lineIndPlot),sizePlot(1) ,sizePlot(2));

% Print metric results
disp('Metric results:')
disp(plotArr)

% Plot metric results
figure
imagesc( unique(xPoints), unique(yPoints), log10(plotArr));
title(['(Scrambled Correlation) log_{10} of: ', message])
colormap(gray);
colorbar
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]')

caxis([2 3]) % Normalize color display (expand as needed)

% Debug: show last area based filter
% figure
% toPlotOrig = imgDataObj{:}.get('CData');
% toPlotOrig(jCluttLinInd) = 0;
% imagesc(toPlotOrig);
% colormap(gray);
% caxis([-55 0]);
