% angleVal = line angle in degrees (90 is straight down)
% maxPoint = point to center line about
% cutLength = how long our cutting lines are, in mm
% XRes = distance in mm between x elements
% YRes = distance in mm between y elements
function normThreshLen = widthOnAngle(cutLength,angleVal, maxPoint, imgData, threshdB, xRes, yRes)

% Determine the number of pixels that we need to grab for this line
xPix = abs(cutLength*cos(angleVal)/(xRes));
yPix = abs(cutLength*cos(angleVal)/(xRes));
LPix = sqrt(xPix^2 + yPix^2);
numPointsHalfLine = round(LPix/2);

% Get data along line
[lineData, endPoint1, endPoint2] = dataFromPointAng(imgData, angleVal, maxPoint(1), maxPoint(2), numPointsHalfLine);

% Find where we first hit -45 dB moving inwards along line data
brightInd = (lineData > threshdB);
aboveThresh = lineData(find(brightInd,1,'first'):find(brightInd,1,'last'));

% plot(aboveThresh);

% Get width of area above -45 dB (in pixels)
aboveThreshPix = numel(aboveThresh);

% Find distance between the endpoints
axDist = abs(endPoint2(1) - endPoint1(1));
latDist = abs(endPoint1(2) - endPoint2(2));

% Get total length of line, in mm
normTotLength = sqrt((axDist*yRes)^2 + (latDist*xRes)^2);

% Calculate total length of above threshold componenet, in mm
% Assume thresLen/totLen = thresPix/totPix => threshLen = totLen * (thresPix/totPix)
normThreshLen = normTotLength*(aboveThreshPix/numel(lineData));
