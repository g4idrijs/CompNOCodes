% Provide angle in degrees 
% indAxStart = starting axial index
% indLatStart = starting lateral index
% Grabs 2*numPointsHalfLine along line
function [lineData, endPoint1, endPoint2] = dataFromPointAng(imgData, angleVal, indAxStart, indLatStart, numPointsHalfLine)

% Set angle of line and number of elements in line
angleVal = angleVal*pi/180; % radians (0 degrees is straight down)

% How far we move with each chunk of length in the line
latIncRate = cos(angleVal);
axIncRate = sin(angleVal);

% Get end point of line
endPointLat1 = indLatStart + numPointsHalfLine*latIncRate;
endPointAx1 = indAxStart + numPointsHalfLine*axIncRate;

% Get another endpoint of the line
endPointLat2 = indLatStart - numPointsHalfLine*latIncRate;
endPointAx2 = indAxStart -numPointsHalfLine*axIncRate;

% Get data between endpoints
rowIndStartEnd = round([endPointAx1, endPointAx2]);
colIndStartEnd = round([endPointLat1, endPointLat2]);
lineData = improfile(imgData,colIndStartEnd, rowIndStartEnd); % Returns not a number if lands outside picture

endPoint1 = [endPointAx1, endPointLat1];
endPoint2 =  [endPointAx2, endPointLat2];
