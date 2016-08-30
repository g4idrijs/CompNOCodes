function A = getSolveTotals( solvePairs,lastIndSolvEqn,N,eqnCellArray,valList)
% Form the matrix, the "A" in Ax = b
% This looks at the rows that DO contain the solvePairs
A = [];

    % Each solvePair contributes two columns to the matrix
    for i = 1:size(solvePairs,1)
        currSolvePair = solvePairs(i,:);        
          
        AcolPair = [];
        
        % Each equation contributes one row of the matrix
        for eqnInd = (lastIndSolvEqn+1):(N-1)
            currEq = eqnCellArray{eqnInd};

            % Find rows that contain the solve pairs
            [rowInclSolveVar,~] = find(ismember(currEq,currSolvePair));

            % Get indices of elements in those rows that are not the solve pairs
            sumInd = setdiff(unique(currEq(rowInclSolveVar,:)), solvePairs);
            
            % We have a problem here because sometimes the solve pairs
            % overlap with eachother!
            
            % Sum over the elements with these indices
            hSum = sum(valList(1,sumInd));
            hpSum = sum(valList(2,sumInd));
            
            AcolPair = [AcolPair; [hSum hpSum]];
        end        
       
        A = [A AcolPair];       
    end

end

