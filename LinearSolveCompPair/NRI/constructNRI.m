function NRICode = constructNRI(intervals)

% Insert zeros to make into an NRI form
nonzeroSpots = [1 1+cumsum(intervals+1)];

% Finish making the vector to shift relative to itself
NRICode = zeros(1,5); 
NRICode(nonzeroSpots) = 1;