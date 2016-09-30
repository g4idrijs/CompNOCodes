% Generate a set of codes given the code pairs that should be neighbors
% (be used in adjacent focus zones)

% pair1 = a complementary code pair
% pair2 = a complementary code pair
% pair1 and pair2 should have good crosscorrelation properties

% numFocX = number of focal zones in x
function codeSet =  codeFromNeigh(pair1, pair2, numFocX, numFocZ)

% Determine code length
if(size(pair1,2) == size(pair2,2))
    N = size(pair1,2);
else
    error('Neighbor code pairs are of different lengths.')    
end

% Stores one code per row
codeSet = zeros(numFocX*2,N);

% We alternate codes in x, and use the same code numFocZ times in z
codeCounter = 1; % Which code we are currently saving
for currXFoc = 1:(numFocX)
    for currZFoc = 1:numFocZ              
        if(mod(currXFoc,2) == 1) 
            codeSet(codeCounter,:) = pair1(1,:);
            codeSet(codeCounter+1,:) = pair1(2,:);
        elseif(mod(currXFoc,2) == 0) 
            codeSet(codeCounter,:) = pair2(1,:);
            codeSet(codeCounter+1,:) = pair2(2,:);
        end
        codeCounter = codeCounter + 2;
    end
end


end


