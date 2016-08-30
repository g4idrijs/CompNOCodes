function [ valList, alreadySet ] = setThreeVarAndSolve( eqnCellArray, N, lastIndSolvEqn)

% Keep track of which variables have been set already
alreadySet = zeros(2,N); % First row = h, second row = h'
% Ex. if h'3 has been set, then set (2,3) = 1.

% Represents value of unknown variable
syms x

for eqn = 1:lastIndSolvEqn    
    currEq = eqnCellArray{eqn};

    % Determine variables to set
    % (h=1 || h' = 2, index of variable)
    % Go through all variables, and only include 3 that aren't set yet
    numSetCurr = 0;
    varToSet = [];
    currVars = unique(currEq);
    for i = 1:2 % Go through h and h' values
        for j = 1:length(currVars) % Go through variables in currEq            
            currVar = currVars(j);
            if(alreadySet(i,currVar) == 0)                
                varToSet = [varToSet; [i,currVar]];
                numSetCurr = numSetCurr + 1;
            end
            if (numSetCurr >= 4)
               break; 
            end
        end
        if (numSetCurr >= 4)
               break; 
        end
    end

    % Set appropriate variables randomly
    % We set all but one randomly, and will solve for the last
    for i=1:(size(varToSet,1)-1)
        currInd = varToSet(i,:);
        valList(currInd(1), currInd(2)) = randn(1);  
        
        alreadySet(currInd(1),currInd(2)) = 1;
    end

    % Determine variable to solve for
    lastVarInd = varToSet(end,:);
    varToSolve = [lastVarInd(1) lastVarInd(2)];

    % Solve for unknown variable
    
    % Go through equation and add up all terms for h
    % h1*h5 is what we want
    % If value is not known, substitute with an "x"
    
    % For each pair in the equation...
        % Look up the values of each variable in the pair
        % If one value is the variable value, then sub in "x"
    % This computes the equation we want to set to zero
    runningTotal = 0;
    for i = 1:size(currEq,1)
        currPair = currEq(i,:);
        for j = 1:2  % j = 1 means work with h, j = 2 means work with h'
            if (alreadySet(j,currPair(1)) == 1)
                firstPairElement = valList(j,currPair(1));
            else
                firstPairElement = x;
            end
            if (alreadySet(j,currPair(2)) == 1)
                secPairElement = valList(j,currPair(2));
            else
                secPairElement = x;
            end
            currPairProd = firstPairElement*secPairElement;            
            
            runningTotal = runningTotal + currPairProd;
          
        end
        
    end
    
    % Set the sum of autocorrelations to zero, and solve
    eqn1 = runningTotal == 0;
    sol = double(solve([eqn1],[x]));

    % Store found value in value list
    valList(varToSolve(1), varToSolve(2)) = sol;
    alreadySet(varToSolve(1),varToSolve(2)) = 1;       
    
end

end

