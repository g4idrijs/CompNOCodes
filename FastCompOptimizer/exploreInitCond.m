clear

 % Where to start loading pairs (we load 20 codes each time - 10 pairs)
loadStarts = 1:10:65;

for i = 1:numel(loadStarts)
    startLoad = loadStarts(i); % Where we are currently loading pairs from    

    % Using the startLoad specified seed, generate and save complementary
    % pairs
    nonLinearGenPairsConstrFastCompVer
end