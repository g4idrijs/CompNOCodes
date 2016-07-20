function A = getSolveTotals( solvePairs,lastIndSolvEqn,N,eqnCellArray,valList)
% Form the matrix, the "A" in Ax = b
% This looks at the rows that DO contain the solvePairs
A = [];
% Start with first eqn for now   
for eqnInd = (lastIndSolvEqn+1):(N-1)
    currEq = eqnCellArray{eqnInd};

    % Find rows that contain the solve pairs
    [rowInclSolveVar,~] = find(currEq == solvePairs);

    % Get indices of elements in those rows that are not the solve pairs
    sumInd = setdiff(unique(currEq(rowInclSolveVar,:)), solvePairs);

    % Sum over the elements with these indices
    hSum = sum(valList(1,sumInd));
    hpSum = sum(valList(2,sumInd));

    % These two sums form another row of the A matrix
    Arow = [hSum, hpSum];
    A = [A; Arow];
end


end

