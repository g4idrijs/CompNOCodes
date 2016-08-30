function [  intVect] = createIntVect( lenVect, seed, issAllow, badInt)
% Start with the interval vector "seed", creates 
% a vector of interval lengths intVect that has an upper bound
% on tag issues of issAllow and length lenVect.

% badInt is dissallowed interval lengths

% The codes in the paper uses ubIssues = 0

% Strategy:
% Start as small as possible, and see if that interval can be 
% added without exceeding ubIssues. If so, add it.
% Continue until length is as long as desired.

% How many tag issues we can add still?
issuesCanUse = issAllow - ubIssues(seed); 

intVect = seed;

% Keep adding intervals until we have enough intervals
while(length(intVect) < lenVect)    
    haveNext = 0; % Do we have the next interval yet?
    i = 1; % The proposed next interval length    
    while(haveNext == 0)
        if (any(i == badInt) == 0) % Avoid excluded intervals
            % Check to see if we can handle the increase in issues
            numIssuesIncr = incUbIssues(intVect, i);
            if (numIssuesIncr <= issuesCanUse)
                intVect = [intVect, i];
                issuesCanUse = issuesCanUse - numIssuesIncr;
                haveNext = 1;
            end
        end
        i = i + 1;
    end        
end


end

