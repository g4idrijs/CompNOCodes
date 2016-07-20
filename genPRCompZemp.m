% Use Dr. Zemp's method to generate complementary pairs
% We randomly pick some elements, and then solve for the rest

%% Generating equations and storing as index pairs
% Ex. h1*h4 + h2*h5 = -(h1'*h4' + h2'*h5')
% is stored as [1 4; 2 5]

% Since each equation is a 2D matrix, we store sets of equations
% as a cell, where each cell is a 2D matrix
N = 5; % The length of the codes
eqnCellArray = cell(1,N-1);

% The j_th equation is given by: (for j = 1..N-1)
% sum(h_i*h_{N+i-j}) for i = 1..j = sum(h'_i*h'_{N+i-j}) for i = 1..j
for j = 1:(N-1)
    currEqn = [];    
    for i = 1:j
        % Current term - an index pair (ex. h1*h4)
        currPair = [i N+i-j];      
        currEqn = [currEqn; currPair];   
    end
    % Current equation (ex.  h1*h4 + h2*h5 = -(h1'*h4' + h2'*h5'))
    eqnCellArray{j} = currEqn; 
end

printCellArray(eqnCellArray)

%% Determine which equations we can solve individually 
% ..by setting 3 variables of the 4 new variables in them

% We set exactly N+1 variables
% We can solve equations individually when we can set 3 more variables
% So, the question is how many times 3 goes into N+1

% The index of the last equation we solve individually
lastIndSolvEqn = floor((N+1)/3); 

%% Create a list to store variables values
% Store h on first row, h' on second row
% Ex. store h_2 on the second entry of the first row
valList = zeros(2,N);

%% Solve the equations that can be solved individually
% Set three new variables to random values
% and then solve for the remaining variable

% Keep track of which variables have been set already
alreadySet = zeros(2,N); % First row = h, second row = h'
% Ex. if h'3 has been set, then set (2,3) = 1.

for eqn = 1:lastIndSolvEqn    
    currEq = eqnCellArray{eqn}

    % Determine variables to set
    % (h=1 || h' = 2, index of variable)
    % Go through all variables, and only include 3 that aren't set yet
    numSetCurr = 0;
    varToSet = [];
    for i = 1:2 % Go through h and h' values
        for j = unique(currEq) % Go through variables in currEq            
            
            if(alreadySet(i,j) == 0)                
                varToSet = [varToSet; [i,j]];
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

    valList
    
    % Determine variable to solve for
    lastVarInd = varToSet(end,:);
    varToSolve = [lastVarInd(1) lastVarInd(2)]

    % Solve for unknown variable
    syms x
    % Go through equation and add up all terms for h
    % h1*h5 is what we want
    % If value is not known, substitute with an "x"
    
    % For each pair in the equation...
        % Look up the values of each variable in the pair
        % If one value is the variable value, then sub in "x"
    
    
    eqn1 = valList(1,1)*valList(1,5) == -valList(2,1)*x;
    sol = double(solve([eqn1],[x]));

    % Store found value in value list
    valList(varToSolve(1), varToSolve(2)) = sol;
    alreadySet(varToSolve(1),varToSolve(2)) = 1;
    
    alreadySet

    valList
    xcorr(valList(1,:))  + xcorr(valList(2,:))
    
    a=3;
end























