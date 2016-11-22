% Plots width down to threshold values at all angles for 
% specified spot in phantom

clear

%% High level parameters
% Plot width at different cross sections for points of interest
% Points of interest
pointsIntList = [0 0 45e-3]; %[0 0 30e-3; 0 0 40e-3; 0 0 50e-3];
numPointsInt = size(pointsIntList,1);

% Threshold for determining width at angles
threshdB = -45; % dB

%% Load data
% Get data from loaded figure
h = gcf;
hChildren = h.Children;  
dataObjs = get(hChildren, 'Children'); %handles to low-level graphics objects in axes
imgDataObj = dataObjs(2); % Not very robust! (Would like to call by name)

% Load figure data
imgData = imgDataObj{:}.get('CData');

% Get axis limits from original graph
currA = gca;
sourceYLim = currA.get('YLim');
sourceXLim = currA.get('XLim');

% Grab old title, in case it's wanted for output
titleVal=get(currA,'Title');
titleVal=get(titleVal,'String');

% Get y axis, x axis values
yAxis = linspace(sourceYLim(1),sourceYLim(2),size(imgData,1))/1000; % m
xAxis = linspace(sourceXLim(1),sourceXLim(2),size(imgData,2))/1000; % m
yRes = (yAxis(2) - yAxis(1))*1000; % mm per pixel in y direction
xRes = (xAxis(2) - xAxis(1))*1000; % mm per pixel in x direction

% Measure widths assuming pixels in each direction are 1 mm
a=3;
% yRes = 1;
% xRes = 1;

%% Loop through points of interest
figure

for phantPointInd = 1:numPointsInt
    %% High level parameters

    % Point of interest (lateral position should be center of phantom)
    phantPointInt = pointsIntList(phantPointInd,:); % (x,y,z) [m]

    %% Find initial line center

    % Get position of point (center for line)
    latStart = phantPointInt(1);
    axialStart = phantPointInt(3);

    % Get indices of start point (center for line)
    [valLatStart,indLatStart] = min(abs(xAxis-latStart));
    [valAxStart,indAxStart] = min(abs(yAxis-axialStart));

    %% Find actual maximum point to center operations

    % Get index of row that is maximum
    [maxVal,maxRowIndOffset] = max(imgData((indAxStart - 50):(indAxStart + 50),indLatStart));
    maxRowInd = (indAxStart - 50) + maxRowIndOffset - 1;

    % Get actual max point
    maxPoint = [maxRowInd, indLatStart];

    %% Take slices and get width that is always above -45 dB
    % (move in from outside until we hit -45 dB)
    % Get data along the line
    
    cutLength = 2; % mm (width of line we use to get through points)    
    numLines = 300;
    angleValRange =  linspace(0,180,numLines); % Angles to determine width at, in degrees. 90 is straight down, and the angles increase clockwise

    aboveThreshWidths = zeros(1,numel(angleValRange));
    angleValInd = 1; % Keeps track of which angle we are storing
    for angleVal = angleValRange        
        aboveThreshWidths(angleValInd) = widthOnAngle(cutLength,angleVal, maxPoint, imgData, threshdB, xRes, yRes);
        angleValInd = angleValInd + 1;
    end

    hold on
    plot(angleValRange,aboveThreshWidths)

end
%% Add a legend
N = 1:numPointsInt;
% legend(cellstr(num2str(N', 'Point %-d')));
title('Width Above Threshold')
xlabel('Angle (degrees)')
ylabel('Width (mm)')
ylim([0 max(aboveThreshWidths)])






