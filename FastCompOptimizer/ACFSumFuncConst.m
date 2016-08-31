% Used when solving autoccorrelation constraint using fmincon

% ceq is a thing we want to set to zero
% c is used for an inequality-based constraint

% We want the ratio of abs(sidelobes)/abs(mainLobe) to be small

% N == code length
% numPairs = number of complementary pairs

function [c,ceq] = ACFSumFuncConst(x,N,numPairs,intervals)

    %c = [maxXcorr(x, intervals)/minMainLobe(x, intervals) - 0.1];
    c = [];
    %x = reshape(x, [5*2, 8]);
    if length(intervals)
        x_int = zeros(size(x, 1), size(x, 2)+sum(intervals(1:size(x, 2)-1)));
        k = 2;
        x_int(:, 1) = x(:, 1);
        for i = 2:size(x, 2)
            insert = [zeros(size(x, 1), intervals(i-1)) x(:, i)];
            x_int(:, k:k+size(insert, 2)-1) = insert;
            k = k+size(insert, 2);
        end
        x = x_int;
    end
        
    % Store the autocorrelation functions of each pair    
    ceq = zeros(size(x, 2)-1,numPairs);
    
    for pair = 1:numPairs
        
        % Get current pair        
        h = x(2*pair-1,:);
        hp = x(2*pair,:);
        
        % Get mainlobe value in sum of autocorrelations
        mainLobeValPair = h*h' + hp*hp';    
        
        % Find sum of ACFs
        totXCorr = xcorr(h) + xcorr(hp);
        
        % Grab the relevant half of the ACF (side lobe)
        % and then divide by main lobe value
        ceq(:,pair) = abs(totXCorr(1:size(x, 2)-1))./mainLobeValPair;
      
    end  

end