function animInfo = visBeamSubModule(xmt,codes,c,fs, plotSingLoc)

% Visualize the beam currently planned for xmt

%% Plot excitation over time at a single point
if(plotSingLoc == 1)
    visLoc = [0 0 40/1000];
    [h,~] = calc_hp(xmt,visLoc);
    
    if(numel(h) > 1000)
       h = h(1:1000); 
    end
    
    figure
    plot(h)
    title(['Beam at [', num2str(visLoc),'] mm During First Transmission']);
    ylabel('Emitted Pressure Field')
    xlabel('Time Index (right is later in time)')
    
    animInfo = {};
    
    % Convolve with the spatial impulse response (get expected echo)
%     [impResp, ~] = calc_h(xmt, visLoc);
%     expecEcho = conv(impResp,h);
%     figure
%     plot(expecEcho)
%     title('Expected Echo')
    
else
    %% Plot excitation over a region
    % Get excitation at a range of positions over time
    xRange = [-10:1:10]/1000; % m
    zRange = [5:0.1:40]/1000; %[35:0.1:45]/1000; % m    
    pointsInterest = zeros(numel(zRange),3);
    pointsInterest(:,3) = zRange';  

    % Determine indices of samples to keep
    [Rmax, Rmin, ~, ~, ~, Smin_cVis, Smax_cVis, ~, no_rf_samples_cVis] =...
    calcSampleTimeRanges(pointsInterest, codes, c, fs);
    Smin_cVis = round(Smin_cVis/2); % Adjust for no echo here
    Smax_cVis = round(Smax_cVis/2);
    no_rf_samples_cVis = round(no_rf_samples_cVis /2);

     % Create matrix to store data (preallocate)
    respMat = zeros(numel(zRange), numel(xRange), no_rf_samples_cVis);

    outputMsg = '';
    fprintf('\n')
    for xRangeInd = 1:numel(xRange)
        % Print current frame
        fprintf(repmat('\b', 1, length(outputMsg)));
        outputMsg = sprintf('Frame %d/%d', xRangeInd, numel(xRange));
        fprintf(outputMsg);            
       for zRangeInd = 1:numel(zRange)
           currX = xRange(xRangeInd);
           currZ = zRange(zRangeInd);
          % Calculate emitted field at particular place over time
           [currResp,start_timeVis] = calc_hp(xmt,[currX 0 currZ]);  

           % Grab needed data and align in time
           currResp = alignRF(currResp,start_timeVis,fs,Smin_cVis,Smax_cVis,no_rf_samples_cVis,1);

           % Store sequence at indices corresponding to (x,z) location
           respMat(zRangeInd,xRangeInd,:) = currResp;
       end
    end

    %% Show the animation
    figure
    numFrames = size(respMat,3);
    for frame = 1000; 
        imagesc([min(xRange), max(xRange)]*1000,[Rmin, Rmax]*1000, respMat(:,:,frame))

        title(['First Code Set. Middle Line. Frame: ', num2str(frame), ' of ', num2str(numFrames)]);
        colorbar;
        caxis(1.0e-12 *[-0.1417,0.1398])

        ylabel('z (mm)')
        xlabel('x (mm)')   

    end   
    
    % Allows us to call up a different frame as desired
    animInfo = {respMat, xRange, Rmin, Rmax};
end

end