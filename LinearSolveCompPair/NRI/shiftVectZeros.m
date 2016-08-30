% For visualizing cross correlations
% when there are zeros in known locations
 
clear

% Construct and display the code
zeroIntervals = [1 2];
v = constructNRI(zeroIntervals);
disp(v);
fprintf('\n');

% Display shifted copies
numShift = length(v) - 1;
for i = 1:numShift   
   disp(v)
   v = [0, v];
end