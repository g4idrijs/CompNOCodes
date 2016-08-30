function ubIssues = ubIssues( intVect )
% Finds maximum number of times a join length is repeated
% is a sequence of lengths
% (an upper bound in the number of tag alignment issues)
% For example, [2,3,2,2] has only 2 tag alignment issues, but this function
% will return 3.

vectJoins = joins(intVect);
[~,f] = mode(vectJoins); % Grabs the number of times the mode shows up

ubIssues = f - 1;

end

