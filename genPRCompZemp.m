% Use Dr. Zemp's method to generate complementary pairs
% We randomly pick some elements, and then solve for the rest

% This code designed with N = -1 mod 3 and N > 2 in mind.
% CURRENTLY ONLY WORKS FOR N = -1 mod 3
% (Should also work with N = 1 mod 3, but not with N = 0 mod 3)

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

%printCellArray(eqnCellArray)

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

[ valList, alreadySet ] = setThreeVarAndSolve( eqnCellArray, N, lastIndSolvEqn);

% Display values found, and sum of autocorrelation so far
valList
xcorr(valList(1,:))  + xcorr(valList(2,:));

%% Print the remaining system of equations
disp('Equations left to solve:')
for i = (lastIndSolvEqn+1):(N-1)
    disp(eqnCellArray{i})
end

% FIRST solve the remaining equations if N = -1 mod 3
% In this case, we don't have to set any more values
% -we have equal number of equations unknowns remaining

%% Determine variables to solve for (DEBUG)
% (h == 1 or h' == 2, location in h or h')
[~,colToSolve] = find(alreadySet == 0);
solvePairs = unique(colToSolve); % Ex. if solvePairs = 3, solve for h3 and h3'

%% Form the column vector 
% ...that is the negative of the sum of the non-solve terms. 
% Is "b" in Ax=b.

nonSolveTotals= getNonSolveTotals(solvePairs,lastIndSolvEqn,N,eqnCellArray,valList);
b = -nonSolveTotals'

%% Form the matrix, the "A" in Ax = b
% This looks at the rows that DO contain the solvePairs
A = getSolveTotals(solvePairs,lastIndSolvEqn,N,eqnCellArray,valList)

%% Solve the matrix equation to get the remaining variables
solveVals = A\b;
valList(:,solvePairs) = solveVals

xcorr(valList(1,:)) + xcorr(valList(2,:))












