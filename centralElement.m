function centralEl = centralElement( v )
% Returns the central element of a vector with an odd number of elements
% For example, for [5,-1,3] it returns -1.
% It errors if an even vector is fed to it.

% Make sure there is a central element
vLen = length(v);
assert(rem(vLen,2) == 1, 'Vector needs to have odd length.')

% Get central element
indCent = (vLen - 1)/2 + 1;
centralEl = v(indCent);
end

