function incUbIssues = incUbIssues(currIntVect, propInt)
% currIntVect = a vector holding the current interval lengths, in order
% (ex. (2,3,4))
% propInt = the proposed next interval
% This function returns the increase in the upper bound number of problems if we add propInt
% This should be a safe upper bound on the increase in the maximum element
% of the cross correlation matrix (some problems will happen at seperate
% times).

% Ex. incUbIssues([2,3],[5]) should return 1, because there is one new
% way we can shift the code with respect to itself to get a dot product
% greater than 1 (two tag bits line up at once).
% Ex.   10100100001
%       _____101001 (lines up in two places) 

% An issue is created if the sum of the lengths of continguous intervals -
% let's call this a "join length" - already present match with the join
% lengths created just because of the added interval.

incUbIssues = ubIssues([currIntVect,propInt]) - ubIssues(currIntVect);

% Note that adding issues doesn't always increase the upper bound on
% issues.
% However, if we keep the upper bound below desired values we are OK.

end

