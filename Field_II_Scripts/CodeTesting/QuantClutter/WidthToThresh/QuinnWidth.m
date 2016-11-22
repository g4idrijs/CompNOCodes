% Created by Quinn
% Edited by David

function rewardVal = QuinnWidth(mouseInp, pointIntVal, plotResultsVal,figureImg)

%% High level paramters

% Get point of interest with mouse
% mouseInp = 0;

plotResults = plotResultsVal; % Plot results
holdPlot = 0; % Add to past results in same figure

% Get point of interest
% m (lateral, axial) (point of interest)
if(mouseInp == 0)
   pointInt = pointIntVal; %[0e-3 40e-3];
else
    pointInt=ginput(1)/1000;
end

% Number of lines to cut through the point
numCutLines = 360;

% Set threshold for determining bright regions, and title
threshdB = -45; % dB
titleStr = sprintf('Walking Aperture: Center, Focused. Thresh: %i', threshdB);

cutLineLen = 2e-3; % m (half-length of cut line: sort of a radius)
maxSearchBoxSideLen = 2e-3; % m (how large of a box about pointInt to search for a maximum)

%% Grab data from loaded figure (set as gcf to load current data)
h = figureImg;
hChildren = h.Children;  
dataObjs = get(hChildren, 'Children'); %handles to low-level graphics objects in axes
try
    imgDataObj = dataObjs(2); % Not very robust! (Would like to call by name)
catch
    error('Please ensure PSF image is in focus.')
end

% Load figure data
imgData = imgDataObj{:}.get('CData');

% Get axis limits from original graph
currA = gca;
sourceYLim = currA.get('YLim');
sourceXLim = currA.get('XLim');

% Get y axis, x axis values (in m)
yAxis = linspace(sourceYLim(1),sourceYLim(2),size(imgData,1))/1000; % m
xAxis = linspace(sourceXLim(1),sourceXLim(2),size(imgData,2))/1000; % m
yMPix = abs(yAxis(2) - yAxis(1)); % size of axial pixels, m
xMPix = abs(xAxis(2) - xAxis(1)); % size of lateral pixels, m

% Make data positive, and remove data below threshold
newImage = imgData - threshdB;
newImage(newImage <= 0) = 0;

sizeZeros = [40 312];   % size of rectangle you're gonna 'zero out' when trying to find maximums
tempImage = newImage;

% finds the maximum, records its place, 'zeros out' all data in a rectangle
% with the centre corresponding to that maximum in order to 'remove dot'
% from image

%% Get coordinates of maximum intensity near point of interest

% Starting point of interest (in m)
latPointInt = pointInt(1);
axPointInt = pointInt(2);

% Start point of interest (indices)
[~,indLatStart]  = min(abs(xAxis-latPointInt));
[~,indAxStart] = min(abs(yAxis-axPointInt));

% Look in a box around this point for a max
latPixBox = round(maxSearchBoxSideLen/xMPix);
axPixBox = round(maxSearchBoxSideLen/yMPix);
boxRowInd = (indAxStart-round(axPixBox/2)):(indAxStart+round(axPixBox/2));
boxColInd = (indLatStart-round(latPixBox/2)):(indLatStart+round(latPixBox/2));
boxData = imgData(boxRowInd,boxColInd);

% Get indices of maximum with respec to box
[peakInt,ind] = max(boxData(:));
[maxRowBox,maxColBox] = ind2sub(size(boxData),ind);

% Finally, get indices of maximum with respect to image
xD = boxRowInd(maxRowBox);
yD = boxColInd(maxColBox);

maxVec = [xD yD];

%% Calculate widths

% Get location of maximum point in meters
maxnumVec = zeros(size(maxVec));
maxnumVec(:,1) =  yAxis(maxVec(:,1));   % converts samples into meters
maxnumVec(:,2) = xAxis(maxVec(:,2));    % converts samples into meters

% Specify radius of cutting lines, in meters
rx = cutLineLen;
ry = rx;

% Specify cutting angles, in radians
theta = linspace(0,pi,numCutLines);

% Will store calculated widths
width = zeros(length(maxVec(:,1)),numCutLines);
lengthvar = width;

% Get data along each cutting line
for lineInd = 1:numCutLines
    % Start point
    xI = rx*sin(theta(lineInd)) + maxnumVec(1,1);
    yI = ry*cos(theta(lineInd)) + maxnumVec(1,2);
    
    % Finish point
    xF = rx*(sin(theta(lineInd) + pi)) + maxnumVec(1,1);
    yF = ry*(cos(theta(lineInd) + pi)) + maxnumVec(1,2);
    
    xpix = [xI xF];
    ypix = [yI yF];
    
    % Get data along line between starting and ending point
    line = improfile(xAxis,yAxis,imgData,ypix,xpix);      
    
    % Find line width
    % If the next dot over is contributing to create intensity
    % over the threshold, that will mess things up.
    wi = find(line > threshdB,1,'first');
    wf = find(line > threshdB,1,'last');
    
    % Normalize width in m
    % (necessary because pixels can be different sizes in different
    % dimensions)
    lengthvar(1,lineInd) = sqrt((xF-xI).^2 + (yF-yI).^2)/length(line); % length of line in m/length in sample.
    width(1,lineInd) = (wf-wi)*lengthvar(1,lineInd); 
end

%% Plot results
% (Note: 90 degrees is straight down)
if(plotResults == 1)
    if(holdPlot == 1)
        figure(5);
        hold on
        plot(theta*(180/pi),width(1,:)*1000); 
        hold off
    else
        figure(5);plot(theta*(180/pi),width(1,:)*1000); 
        ylabel('Width (mm)')
        xlabel('Angle (degrees)')
        title('PSF Width')
        title(titleStr)
    end
end

%% Compare to reference
refPSF = load('refPSF.mat');
refPSF = refPSF.('refPSFWidth');

% Get average diameter of PSF
currPSF = width(1,:)*1000; % mm
blobSize = mean(currPSF); 

% Fraction of global maximum intensity
peakIntVal = 10^(peakInt/20); 

% Reward small blobs with high intensity
rewardVal = peakIntVal/blobSize;



