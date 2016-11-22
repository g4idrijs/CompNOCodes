% Given a set of codes, creates and image and calculates
% the metric at the center point: signal / clutter ratio

% We want to maximize metric

function metric = imageAndGetMetric(codes)

% Create image using codes
addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\CodeTesting\FieldIISim\')
figureImg = codeAndCompSameTime_Function(codes);

% Measure metric of image
addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\CodeTesting\QuantClutter\WidthToThresh\')
mouseInp = 0; pointIntVal = [0e-3 40e-3]; plotResultsVal = 0;
metric = QuinnWidth(mouseInp, pointIntVal, plotResultsVal, figureImg);

% Close invisible figures
close all



