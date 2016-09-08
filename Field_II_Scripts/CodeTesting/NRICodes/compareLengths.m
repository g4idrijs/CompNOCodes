% Run through a variety of code sizes (number of intervals)
% and see how the code lengths compare between:
% -Every interval allowed (or not)
% -2 in sidelobes allowed (vs 1)

seed = [2,3,4];
lenTraditional = zeros(1,100);
lenRelaxed = zeros(1,100);

sizeRange = 4:25;
for size = sizeRange   
    lenTraditional(size-3) = sum(createIntVect(size,seed,0,[1]));
    lenRelaxed(size-3) = sum(createIntVect(size,seed,1,[1]));
end

lenTraditional(lenTraditional == 0) = [];
lenRelaxed(lenRelaxed == 0) = [];

figure
plot(sizeRange,lenTraditional,'o')
hold on
plot(sizeRange,lenRelaxed,'o')
legend('Traditional (bound of 1)','Relaxed (bound of 2)','Location','northwest')
xlabel('Number of intervals between non-zero bits')
ylabel('Sum of interval lengths')
title('Corrected: Code Length Reduction by Relaxing Side Lobe Constraints')