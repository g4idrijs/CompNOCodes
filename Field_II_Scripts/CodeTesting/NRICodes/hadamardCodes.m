function hadMat = hadamardCodes( intVect )
% Returns a matrix, where each row is a code
% The tag spacing is given by the vector intVect
% The coefficients are given by Hadamard matrix coefficients\

% Where we are to place the tags
tagMat = tagMatrix(intVect);

% Create coefficients using Hadamard matrix
coeffMat = hadamard(size(tagMat,1));

% Find where we have tags
tagLoc = find(tagMat);

% Go trhough the tag locations and put in the Hadamard value
hadMat = zeros(size(tagMat));
for i = 1:length(tagLoc)
    currTagLoc = tagLoc(i);
    hadMat(currTagLoc) =coeffMat(i);
end


end

