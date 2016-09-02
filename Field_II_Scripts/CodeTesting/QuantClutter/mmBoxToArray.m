% Specify a box in mm and get the corresponding array data
% imgDataObj = the image data object (from figure)
% boxCenter = center of box (x,y) in mm
% boxWidth = width of box (x) in mm
% boxHeight = height of box (x) in mm
function [arrData,boxRowInd,boxColInd] = mmBoxToArray(imgDataObj, boxCenter, boxWidth, boxHeight)

% Convert box center to a position in array elements (row, column)
% Also convert box dimension into elements
[arrPos,mmPerElsX,mmPerElsY]  = mmToArray(imgDataObj, boxCenter);
widthLines = boxWidth / mmPerElsX;
heightLines = boxHeight / mmPerElsY;

% Get the box data
boxRowInd = arrPos(1) +(-floor(heightLines/2):floor(heightLines/2)); 
boxColInd = arrPos(2) +(-floor(widthLines/2):floor(widthLines/2));
imgData = imgDataObj{:}.get('CData');
arrData = 10.^(imgData(boxRowInd, boxColInd)./20);

    
