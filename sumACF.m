function sumACF = sumACF( codes )
% Returns the sum of the autocorrelation functions
% of the rows of the matrix "codes"


numCodes = size(codes,1);
sumACF = xcorr(codes(1,:));
for i = 2:numCodes
    sumACF = sumACF + xcorr(codes(i,:));
end

end

