% Print the number of variables present in the
% j-th autocorrelation equation for codes of length N.

resultMat = [];

NMax = 6;
Nrange = 2:NMax;
for N = Nrange
    nRow = [];
    for j = 1:NMax
        if (j <= N-1)
            if (j < (N+1)/2)
                nRow = [nRow, 4*j];  
            else
                if (mod(N,2) == 1)
                    nRow = [nRow, 4*floor(N/2)+2]; % Odd case
                else
                    nRow = [nRow, 4*floor(N/2)]; % Even case
                end
            end
        else
            nRow = [nRow, 0];
        end
    end   
    resultMat = [resultMat; nRow];
end
resultMat = [(Nrange)', resultMat];
resultMat