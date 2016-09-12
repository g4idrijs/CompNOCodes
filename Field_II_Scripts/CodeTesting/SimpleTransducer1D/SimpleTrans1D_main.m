%% Simulate output of transducer array (1D)
% Set impulse response, excitation: see output of transducer
clear

%% Specify transducer impulse response and excitation
% We specify the values of two sequences at:
% t = 0, T_s, 2*T_s, 3*T_s where T_s is the period with which the
% transducer can change what it fires
% (we are limited to causual output in this framework)

% Impulse response of transducer on transmission (starting at t = 0)

lenImpResp = 20; % Period, in units of T_s (of impulse response)
hTx = sin(2*pi/(lenImpResp)*(0:lenImpResp)).*hanning(lenImpResp + 1)';

plotImp = 1;
if(plotImp == 1)
    figure
    plot(0:(numel(hTx) - 1), hTx,'o')
    title('Impulse Response')
end

% Excitation (starting at t = 0)%
%exc = [1 zeros(1,lenImpResp-1) 2];  % Space out code
%exc = [1 ones(1,lenImpResp-1) 2 2*ones(1,lenImpResp-1)];    % Repeat code
exc = 1;

%% Calculate transducer output
txOut = conv(hTx,exc);

%% Display transducer output
figure
subplot(1,4,1)
plot((1:numel(txOut)) - 1,txOut)

xlabel('Time')
ylabel('Magnitude')

% Show excitation as text below the graph
str = '';
for i = 1:length(exc)
    if (i == 1)
        str  = strcat(str, sprintf('\n Excitation: [%1.1f',exc(i)));  
    else
        str  = strcat(str, sprintf(',%1.1f',exc(i)));
    end
end
endstr = strcat(str, ']');
title(strcat(sprintf('Transducer Output \n Impulse Response: Single Cycle Sine Wave'),endstr))

%% Calculate signal received after going through impulse response again
hRc = hTx; % Specify impulse response on receive
rcvIn = conv(hRc, txOut);

% Plot received electrical signal
subplot(1,4,2)
plot((1:numel(rcvIn)) - 1,rcvIn)
title(sprintf('Received at Transducer'))
xlabel('Time')
ylabel('Magnitude')

%% Calculate decoded symbol
decodedSig = xcorr(exc,rcvIn);
subplot(1,4,3)
plot(decodedSig)
title(sprintf('Decoded Signal'))
xlabel('Time')
ylabel('Magnitude')

%% Calculate decoded symbol
absDec = abs(decodedSig);
subplot(1,4,4)
plot(absDec)
title(sprintf('Absolute Value Decoded Signal'))
xlabel('Time')
ylabel('Magnitude')