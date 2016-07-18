% Following "Complementary Sets of Sequences"
% by Tseng and Liu

% Given a complementary orthogonal matrix, we generate a larger one
% We do this in two ways

%% We begin by constructing an complementary orthogonal matrix

% We start with a complementary pair
% (each row is a code)
% A =[2  1  2 -1
%     2 -1 -1  1];

A = [-1 -1  0
    0  1 -1];

% We generate a pair that is a mate to A
% That is, the pairwise cross correlation of A and Amate sums to zero
Amate = [fliplr(A(2:2:end,:))
    -fliplr(A(1:2:end,:))];

% The column chunks correspond to each complementary pair
% ex. if A =[0 -1 -1  1; 0 -1  0 -1]; the chunks are 4 numbers wide.

% Each column chunk is self-complementary (forms a complementary sequence).
% Each column chunk is orthogonal to ecah column chunk 
% (sum of pairwise cross correlations is zero)
   
% Our starting complementary orthogonal matrix
startCO = [A Amate]

% Debug - makes sure chunks have desired properties
dispCheck = 0;
if(dispCheck == 1)
    chunkWidth = size(A,2);
    chunkCol1 = startCO(:,1:chunkWidth);
    chunkCol2 = startCO(:,(chunkWidth+1):2*chunkWidth);
    
    % Check self-complementarity of chunks
    acfChunk1 = sumACF(chunkCol1)
    acfChunk2 = sumACF(chunkCol2)
    
    % Check cross correlation of chunks
    ccf1And2 = sumCCF(chunkCol1,chunkCol2)
end

%% First method of extension
% Number of codes per focus location goes up by power of 2 each time
% Number of focus locations also goes up by power of 2 each time
% I think the patent used this method for generation

% This involves interleaving startCO with itself
% Ex. [1 1   interleaved with -itself is: [1 -1 1 -1
%      1 -1]                               1 -1 -1 1]

% If x = startCO and # means interleave, we need: 
% [x#x (-x)#x
% (-x)#x  x#x]

% We construct the different different interleavings separately

% Constructing startCO#startCO:
% This simply repeats each column twice
% https://www.mathworks.com/matlabcentral/answers/80607-repeat-every-element-in-matrix
%startCO_startCO = reshape(repmat(reshape(startCO',[],1),1,2)',[],size(startCO,1))';
startCO_startCO = kron(startCO,[1 1]);

% Constructing (-startCO)#startCO:
% temp = repmat(reshape(startCO',[],1),1,2)';
% temp(1,:) = -temp(1,:);
% minStartCO_startCO = reshape(temp,[],size(startCO,1))';
minStartCO_startCO = kron(startCO,[-1 1]);

% We can construct the final larger matrix 
% (which is still complementary and orthogonal, with chunks twice as wide)
finalCO = [ startCO_startCO         minStartCO_startCO
            minStartCO_startCO      startCO_startCO];

% Check the matrix has the desired properties
dispChunks = 1;
firstIsGood = checkCO( startCO, finalCO, dispChunks)

%% Second method of extension - doesn't work
% In this method we extend startCO in another way (p. 650)
% -This method also doubles the transmit events each time used

finalCO2 = [kron(startCO,[1 1]) -kron(startCO,[1 1])
            -kron(startCO,[1 1]) kron(startCO,[1 1])];
dispChunks = 0;
secIsGood = checkCO( startCO, finalCO2, dispChunks);
