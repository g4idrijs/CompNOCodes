% Following "Complementary Sets of Sequences"
% by Tseng and Liu

% Given a complementary set of sequences, find another complementary
% set of sequences that is a mate to the first one.
% --This means that the sum of the pairwise cross-correlation between those
% two is zero.

% Note that this only works for even length sequences
% It seems to work at least for <0,-1,1> codes.

% Method on p. 648
% Start with (A1,A2...,Ap)
% Mate is: ([~A2,-~A1],[~A4,-~A3],...])
% where ~ denotes reversal, and - denotes negation

% A complementary set of sequences:
% (Each row is one code)
A =[0 -1 -1  1
   0 -1  0 -1];
acfA = sumACF(A);

mate = [fliplr(A(2:2:end,:))
    -fliplr(A(1:2:end,:))]

% Check ACF and CCF propertis are good
mateCompl = areCompl(mate)
mateCCF = max(abs(sumCCF(A, mate))) < 1e-5

% In general, our sequence will be one column of a matrix, and we generate
% another column by using this algorithm

%% For a complementary pair, there are only two mates
% The second one is the negative of the first:
secMate = -mate

% Check ACF and CCF properties are good with respect to original
mate2Compl = areCompl(secMate)
mate2CCF = max(abs(sumCCF(A, secMate))) < 1e-5



