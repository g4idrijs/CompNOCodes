function sumCCF = sumCCF(codes1, codes2)
% Returns the sum of the pairwise cross corrleation of codes1 and codes2

assert(size(codes1,1) == size(codes2,1), 'Number of codes in codes1 and codes2 don''t match');
numCodes = size(codes1,1);

sumCCF = xcorr(codes1(1,:), codes2(1,:));
for i = 2:numCodes
    sumCCF = sumCCF + xcorr(codes1(i,:), codes2(i,:));
end

end

