% For visualizing cross correlations

v = 1:6;

numShift = length(v) - 1;

disp(v);
for i = 1:numShift
   v = [zeros(1,1) v];
   disp(v)
end