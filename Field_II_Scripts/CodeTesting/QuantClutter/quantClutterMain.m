% clear

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

% Display areas used for analysing image as boxes on figure or no
showBox = 1;

% Use scripted selection of signal box and noise box or no
scriptBox = 0;

% User draws only signal box, and clutter box extends vertically
vertClutter = 1;

% Use just center line
justCentLine = 1;

% Get data from loaded figure
h = gcf;
hChildren = h.Children;  
dataObjs = get(hChildren, 'Children'); %handles to low-level graphics objects in axes
imgDataObj = dataObjs(2); % Not very robust! (Would like to call by name)


if(justCentLine == 1)
    keepCompr = 1; % Remove log compressin
    
    % Load figure data
    imgData = imgDataObj{:}.get('CData');
    if(keepCompr == 1)
        arrData = imgData;
    else
        arrData = 10.^(imgData./20); % Remove compression
    end     
    
    % Get axis limits from original graph
    currA = gca;
    sourceYLim = currA.get('YLim');
    
    % Get central column
    midColIndex = round(size(arrData,2)/2); 
    centLineData = arrData(:,midColIndex);
    
    figure
    
    plot(linspace(sourceYLim(1),sourceYLim(2),numel(centLineData)),centLineData);
    
    codeLen = 3;
    numRep = 1;
    title(sprintf('1 Focal Zone. Central Vertical Line of Point Spread Function \n Code length: %d. Times repeated: %d.', codeLen, numRep))
    xlabel('Axial Distance (mm)')
    ylabel('Normalized Intensity (dB)')
    xlim(sourceYLim);
   
else

%% List points to investigate
analysisPoints  = [ 0 0 40          
%            0 0 30; 0 0 35; 0 0 40; 0 0 45; 0 0 50;           
%            -6 0 30; -6 0 35; -6 0 40; -6 0 45; -6 0 50;
%            6 0 30; 6 0 35; 6 0 40; 6 0 45; 6 0 50;
           ]; % (x,y,z) (mm)
analysisPoints = [analysisPoints(:,1) analysisPoints(:,3)];

    for pointCount = 1:size(analysisPoints,1)
               
        % Use script specified signal and clutter box dimensions
        if(scriptBox == 1)
            %% Set region of signal and of clutters 
            center = analysisPoints(pointCount,:); % (x,z) in mm
            
            % Signal box
            sigWidth = 0.3; % mm
            sigHeight = 0.4; % mm
            sigCentMm = center + [0 sigHeight/2]; % (x,z) in mm

            % Clutter box
            clutterWidth = 4; % mm
            clutterHeight = 3; % mm
            clutterCentMm = center + [0 sigHeight/2]; % (x,z) in mm
            

        else
            % User selects signal and clutter boxes
            sigData = getrect; % [xmin ymin width height] (x axis ->, y axis goes up page)
            sigWidth = sigData(3); % mm
            sigHeight = sigData(4); % mm            
            sigCentMm = [analysisPoints(pointCount,1), sigData(2) + sigHeight/2];
            
            if (vertClutter == 1)
                % Extend signal box to clutter box, vertically
                extFactor = 4; % Factor of vertical extension
                clutterWidth = sigWidth; % mm
                clutterHeight = extFactor*sigHeight; % mm            
                clutterCentMm = sigCentMm;
            else            
                clutterData = getrect; % [xmin ymin width height] (x axis ->, y axis goes up page)
                clutterWidth = clutterData(3); % mm
                clutterHeight = clutterData(4); % mm            
                clutterCentMm = [clutterData(1) + clutterWidth/2, clutterData(2) + clutterHeight/2];
            end
        end
               
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
        
        for metrInd = 1:2
            % Will store metric results
            % One result for each point
            metricArr = zeros(size(analysisPoints,1),1);
            if(metrInd == 1)
                [metric, message ] = maxToMax( arrDataSig, justClutterData); 
            elseif(metrInd == 2)
                [metric, message ] = maxToAvg( arrDataSig, justClutterData);
            end   
            metricArr(pointCount) = metric;
            
            % Print metric results
            disp(message)
            disp(metric)
        end

    end

%     % Visualize the metrics handed out as an array
%     % We want to plot this in positions corresponding to the analysis point
%     % positions
%     yPointsChunk = 5; % mm
%     xPointsChunk = 6; % mm
% 
%     xPoints = analysisPoints(:,1);
%     yPoints = analysisPoints(:,2);
% 
%     xChunkInd = xPoints./xPointsChunk;
%     xChunkInd = xChunkInd - min(xChunkInd) + 1;
% 
%     yChunkInd = yPoints./yPointsChunk;
%     yChunkInd = yChunkInd - min(yChunkInd) + 1;
% 
%     sizePlot = [numel(unique(yPoints)), numel(unique(xPoints))];
%     lineIndPlot = sub2ind(sizePlot,yChunkInd,xChunkInd);
% 
%     plotArr = reshape(metricArr(lineIndPlot),sizePlot(1) ,sizePlot(2));
% 

%     % Plot metric results
%     if(numel(plotArr) > 1)
%         figure
%         imagesc( unique(xPoints), unique(yPoints), log10(plotArr));
%         title(['(Scrambled Correlation) log_{10} of: ', message])
%         colormap(gray);
%         colorbar
%         xlabel('Lateral distance [mm]');
%         ylabel('Axial distance [mm]')
% 
%         caxis([2 3]) % Normalize color display (expand as needed)
%     end
% 
% end

% Show last analysis regions as boxes on graph
if(showBox == 1)
    f1=figure;

    toPlotOrig = imgDataObj{:}.get('CData');
    toPlotOrig(jCluttLinInd) = 0;

    imagesc(toPlotOrig);
    caxis([-55 0]);
end
end