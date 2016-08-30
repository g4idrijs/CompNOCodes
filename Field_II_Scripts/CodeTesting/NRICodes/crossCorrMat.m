function ccMat = crossCorrMat( codeMat )
% Returns the cross correlation matrix of the codes in codeMat
% Each row of codeMat is another code
% crossCorrMat(i,j) is the maximum cross correlation value of
% codes i and j.

% Go through each pair and find maximum value in cross correlation
ccMat = zeros(size(codeMat,1));
for i = 1:size(ccMat,1)
    for j = 1:i                    
        ccMat(i,j) = maxCrossCor(codeMat(i,:), codeMat(j,:)); 
        ccMat(j,i) = ccMat(i,j); % By symmetry
    end
end


end

