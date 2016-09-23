function showFrame(frameList, animInfo)    
% Unpack animation information
% animInfo = {respMat, xRange, Rmin, Rmax};
respMat = animInfo{1};
xRange = animInfo{2};
Rmin = animInfo{3};
Rmax = animInfo{4};

%% Show the animation at a list of provided frames
figure
numFrames = size(respMat,3);
for frame = frameList; 
    imagesc([min(xRange), max(xRange)]*1000,[Rmin, Rmax]*1000, respMat(:,:,frame))

    title(['First Code Set. Middle Line. Frame: ', num2str(frame), ' of ', num2str(numFrames)]);
    colorbar;
    caxis(1.0e-12 *[-0.1417,0.1398])

    ylabel('z (mm)')
    xlabel('x (mm)')   
    
    pause(eps);
end  