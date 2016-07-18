function  isGood = checkCO( startCO, finalCO, dispChunks)

isGood = 1;

% Construct the column chunks
colChunkWidth = size(startCO,2);

for i = 1:colChunkWidth:size(finalCO,2)
    currChunk = finalCO(:,i:(i+colChunkWidth-1));
    if (dispChunks == 1)
        disp(currChunk);
    end
    
    % Outputs a 1 if the current chunk is self-complementary
    isSelfCompl = areCompl(currChunk);
    
    if(isSelfCompl ~= 1)
       isGood = 0;
    end
    
    
    % Output a 1 if currChunk and nextChunk have pairwise cross
    % correlations that sum to zero
    for j = (i+colChunkWidth):colChunkWidth:size(finalCO,2) 
        nextChunk = finalCO(:,j:(j+colChunkWidth-1));
        goodCC = max(abs(sumCCF(currChunk, nextChunk))) < 1e-5; 
        
        if(goodCC ~= 1)
            isGood = 0;       
        end
    end
end

end

