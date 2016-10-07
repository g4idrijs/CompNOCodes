% Returns maximum of the following, over all pairs i:
% max(xcorr(pairs_i, pairs_i^c))

% x = a matrix where each two rows forms a complementary pair
function maxSelfCC = maxSelfCCFun(x)

    % Store largest cross correlation between pairs
    maxSelfCC = -1;

    numRows = size(x,1);
    if(mod(numRows,1) == 0) % Make sure we have complete code pairs
        numCodes = numRows / 2;
        for currCode = 1:numCodes
           code1 = x(currCode*2-1,:);
           code2 = x(currCode*2,:);
           currMaxXcorr = max(abs(xcorr(code1,code2)));
           
           % Update the largest cross correlation found so far
           if(currMaxXcorr > maxSelfCC)
              maxSelfCC = currMaxXcorr; 
           end
        end
    else
       error('Expected even number of codes in matrix x.') 
    end


end

