function areCompl = areCompl( codeSet)
% Returns 1 if the codes in codeSet are complementary, and 0 otherwise
% A code set is complementary if the sum of their autocorrelations is 
% in the shape of an ideal (discrete) delta function.

% Each row in codeSet is another code
    sumAutCorr = xcorr(codeSet(1,:),codeSet(1,:));
for i = 2:size(codeSet,1)
    sumAutCorr = sumAutCorr + xcorr(codeSet(i,:),codeSet(i,:));
end

numNonZero = length(sumAutCorr(abs(sumAutCorr)>0.001));

% Check the resulting sum of autocorrelations looks like a delta function
if (centralElement(sumAutCorr) == 0)
    areCompl = 0; % Need to have a peak in the middle 
elseif (numNonZero ~= 1)
    areCompl = 0;
else
    areCompl = 1;
end
    

end

