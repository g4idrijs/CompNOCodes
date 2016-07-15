% Demo new sequence generation from DFT matrix
% Following "OFDM Codes for Peak to Average Power Reduction And Error
% Correction"

% Begin with a complementary pair
A = [0 -1  0 -1  0];
B = [0  0 -1  0  1];

% Multiply by a column of the corresponding DFT matrix
% (Mmatrix that carries out transform)
T = dftmtx(length(A));

colChoice = 2;
ANew = A.*T(:,colChoice)'
BNew = B.*T(:,colChoice)'
xcorr(ANew) + xcorr(BNew)

% This in general creates complex codes