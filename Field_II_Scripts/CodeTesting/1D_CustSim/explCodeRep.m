% Try repeating code elements diffferent fractions of impulse response
% length
len = 600; % Number of nonzero elements in transducer impulse response

ind = 1; % Keep track of subplot position
figure

for custTimes = len*linspace(1/100, 1, 8) % Try a range of number of times to repeat code elements
    
    % Calculate and plot results
    subplot(2,4, ind);    
    oneDimCustMain(round(custTimes), len)  
    
    title(sprintf('NumRepeat/len(ImpResp): %.2f',custTimes/len))

    ind = ind + 1;
end