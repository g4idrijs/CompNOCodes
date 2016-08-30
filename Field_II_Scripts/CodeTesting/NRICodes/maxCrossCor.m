function maxVal = maxCrossCor( v1,v2 )
% Return the maximum interesting value in the cross correlation sequence
% (for keeping cross correlation low)

% If v1 = v2, returns maximum sidelobe value

if (v1 ~= v2)
    maxVal = max(xcorr(v1,v2));
else
    corr = xcorr(v1,v2);
    % Locate and remove central element
    maxInd = (length(corr) - 1)/2 + 1;    
    sideLobes = corr;
    sideLobes(maxInd) = [];
    
    maxVal = max(sideLobes);
end

end

