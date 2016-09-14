h = [1 2 3];
conv(h,h)

% Shift h
shiftBy = 3;
hShift = [zeros(1,shiftBy) h];
conv(hShift, hShift)