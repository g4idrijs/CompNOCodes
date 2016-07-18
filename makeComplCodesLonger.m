% This method takes in a complementary set and an orthogonal matrix, and
% outputs a complementary set, where the individual elements are longer.
% What we really want is MORE codes

% Note that this method should work with complementary sets that have more
% than two elements (we are not restricted to pairs)

% We need an orthogonal matrix (Q*Q^T = k*I) [ex. Hadamard]s
% and we need a complementary set of sequences.

% This provides a template for how to arrange the codes
% Hadamard matrices are orthogonal
% The number of columns is the number of codes here
orthMat = hadamard(4);

% I believe Hadamard matrices are also complementary
% Each row is a code
compSet = hadamard(4);

% Now we use the orthogonal matrix as a guide how to
% multiply and combine the complementary set of sequences.
% Ex. H = [1 1; 1 -1], then we want [A B; A -B]

codeLen = size(compSet,2);
numCodes = size(compSet,1);
assert(numCodes == size(orthMat,2),'orthMat size and number of codes don''t match.')
newCompCodes = zeros(size(orthMat,1), codeLen*numCodes);

% We add the contribution to each code one code at a time
currStart = 1;
lenSec = codeLen*size(orthMat,1); % Size of matrix due to this code
currEnd = lenSec;
for i = 1:numCodes    
    newCompCodes(currStart:currEnd) =  kron(orthMat(:,i),compSet(i,:));   
    currStart = currEnd + 1;
    currEnd = currEnd + lenSec;
end

newCompCodes

% The result should be a new complementary set of sequences (row-wise)
% Let's check the ACF summation property:
numNewCodes = size(newCompCodes,1);
ACFSum = xcorr(newCompCodes(1,:));
for i = 2:numNewCodes
    ACFSum = ACFSum + xcorr(newCompCodes(i,:));
end
ACFSum






