% F = final data to find central width
% (width of central area above dbThresh dB (ex. -45 dB))
% (in units of impulse response length)
function width = oneDWidth(F,len,dbThresh)

% Strategy: Find maximum and work outwards until we hit >-45 dB

% Normalize and conver to dB
normF = (F./max(max(abs(F))));
dBF = 20*log10(abs(normF) + eps);

% Get index of maximum value
[~,maxInd] = max(dBF);

outsideCent = 1; % Keep track of number of elements outside center

% Go to smaller indices until we enter dB threshold
foundEntBig = 0; % Have we moved into high intensities for large indices?
nextTestBig = numel(dBF); % Next index to test - large indices

while(foundEntBig == 0 && nextTestBig >= 1)
    nextValTest = dBF(nextTestBig);
    if(nextValTest > dbThresh)
        foundEntBig = 1;
    else
        foundEntBig = 0;
        outsideCent = outsideCent + 1;
        nextTestBig = nextTestBig - 1;
    end
end

% Go to larger indices until we enter dB threshold
foundEntSmall = 0; % Have we moved into high intensities for small indices?
nextTestSmall = 1; % Next index to test - small indices

while(foundEntSmall == 0 && nextTestSmall <= numel(dBF))
    nextValTest = dBF(nextTestSmall);
    if(nextValTest > dbThresh)
        foundEntSmall = 1;
    else
        foundEntSmall = 0;
        outsideCent = outsideCent + 1;
        nextTestSmall = nextTestSmall + 1;
    end
end

% Calculate width of high intensity center
% Normalize for fineness of simulation (len imp response is constant)
width = (numel(dBF) - outsideCent)/len;