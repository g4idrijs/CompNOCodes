function nonSolveTotals = getNonSolveTotals(solvePairs,lastIndSolvEqn,N,eqnCellArray,valList)

assert(length(solvePairs) == 1, 'Not implemented solvePairs length.')

% Keep track of the non-solve totals (form the "b" in Ax=b)
nonSolveTotals = zeros(1,length(solvePairs)*2);

numEqnDone = 0;
for eqnInd = (lastIndSolvEqn+1):(N-1)
    % Find the terms that contain the elements we are solving for
    currEq = eqnCellArray{eqnInd};
    [rowInclSolveVar,~] = find(currEq == solvePairs);
    
    % Go through each non-solving pair in the current equation
    % and compute the sum of these non-solving pairs
    nonSolvTot = 0;
    noSolveVarRows = setdiff(1:size(currEq,1), rowInclSolveVar);
    for i = 1:length(noSolveVarRows) 

        % Get the current pair indices from the current equation
        currRow = noSolveVarRows(i);    
        currRowInd = currEq(currRow,:);

        % Find the h term (ex. h2*h4)
        firstHVal = valList(1,currRowInd(1));
        secHVal =valList(1,currRowInd(2));
        currHValProd = firstHVal*secHVal;

        % Find the h' term (ex. h'2*h'4)
        firstHPVal = valList(2,currRowInd(1));
        secHPVal =valList(2,currRowInd(2));
        currHPValProd = firstHPVal*secHPVal;

        % Sum the two: h2*h4 + h2'*h4'
        currTotSum = currHValProd + currHPValProd;   
        nonSolvTot = nonSolvTot + currTotSum;     
    end
    numEqnDone = numEqnDone + 1;
    nonSolveTotals(numEqnDone) = nonSolvTot;
end


end

