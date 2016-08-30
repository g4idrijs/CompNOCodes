% Calculate the number of overlaps as we shift a vector relative to itself
% (the number of places where nonzero entires line up in the shifted
% version and in the original unshifted version)

clear

% Construct an NRI code
zeroIntervals = [1 2 3];
v = constructNRI(zeroIntervals); % Is a vector of 1's and 0's
disp(v)

numOverLaps = round(xcorr(v)); % Magnitude = number of overlaps
numOverLaps = fliplr(numOverLaps(1:(numel(v) - 1))); % Grab relevant half

shiftVect = 1:numel(numOverLaps); % How far we shifted to get this overlap

% Display shifts (left column), and resulting overlaps (right column)
disp([shiftVect' numOverLaps'])