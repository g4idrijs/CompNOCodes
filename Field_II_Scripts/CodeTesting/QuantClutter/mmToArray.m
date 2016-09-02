% Convert an image position in mm's to a position in the array
function [arrPos,mmPerElsX,mmPerElsY] = mmToArray(imgDataObj, mmPos)

% Grab image information
imgData = imgDataObj{:}.get('CData');

% Grab axis information
yAxis = imgDataObj{:}.get('YData');
yRange = [yAxis(1) yAxis(end)];

xAxis= imgDataObj{:}.get('XData');
xRange = [xAxis(1) xAxis(end)];

% Find conversion factor between mm and array rows
mmPerElsX = (xRange(2) - xRange(1))/size(imgData,2);
mmPerElsY = (yRange(2) - yRange(1))/size(imgData,1);

% Find position of central region in data array
mmFromLeft = (mmPos(1) - 0) + (0 - xAxis(1));
elsFromLeft = mmFromLeft / mmPerElsX;

mmFromTop = mmPos(2) - yAxis(1);
elsFromTop = mmFromTop / mmPerElsY;

arrPos = round([elsFromTop, elsFromLeft]);

