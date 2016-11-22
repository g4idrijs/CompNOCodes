% Calculate signal to clutter for a set of points

clear

% Points given as (lateral,axial) (m)
numLat = 3;
numAx = 5;
points = [
    -5e-3 30e-3; 0e-3 30e-3; 5e-3 30e-3;
    -5e-3 35e-3; 0e-3 35e-3; 5e-3 35e-3;
    -5e-3 40e-3; 0e-3 40e-3; 5e-3 40e-3;
    -5e-3 45e-3; 0e-3 45e-3; 5e-3 45e-3;
    -5e-3 50e-3; 0e-3 50e-3; 5e-3 50e-3];

rewardVect = zeros(size(points,1),1);


for pointInd = 1:size(points,1)
    disp(sprintf('Point %i of %i',pointInd,size(points,1)));
    
    mouseInp = 0; % Specify input using mouse
    plotResults = 0; % Plot the results
    
    currPoint = points(pointInd,:);
    rewardVect(pointInd) = QuinnWidth(mouseInp, currPoint, plotResults);
end

%%
shapedReward = reshape(rewardVect,[numAx, numLat]);
disp('Reward:')
disp(shapedReward);
