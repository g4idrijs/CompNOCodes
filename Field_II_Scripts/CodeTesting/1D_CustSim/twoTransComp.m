% Return raw decoded data given a system transfer function
% and the side clutter.
% hT == system transfer function
% c == side clutter sequence
% trLen == length of impulse response of transducer
function R = twoTransComp( hT, c, trLen, custTimes, codes)
    % Set codes to use
    % Specify sequence values, starting at n = 0.
    % Unspecified values are assumed to be zero.
    code1 = codes(1,:);
    code2 = codes(2,:);  

    if(custTimes ~= -1)
        % Repeat code elements custom number of times
        numRepeat = custTimes;
    else
        % Repeat code elements (bandwidth matching)
        numRepeat = round(trLen/2);  
    end
    code1 = repelem(code1,numRepeat);
    code2 = repelem(code2,numRepeat);
    
    % Define system function including clutter
    htC = @(inV)conv(hT,inV) + c(length(conv(hT,inV)))';
    
    % Calculate raw decoded data
    R = xcorr(htC(code1),code1) + xcorr(htC(code2),code2);
    
    % Trim zeros (very small values)
    minVal = 1e-10;
    [~,J] = find(R > minVal);
    R = R(min(J):end);
    [~,J] = find(R > minVal);
    R = R(1:max(J));
end

