clear
% Impulse response of transducer
L = 15;
n = 0:L-1;
m = 3;
h = sin(2*pi*m/(L-1)*n);
% h = [0 1 2 0]; 
L = numel(h);

% Self convolution (proportional to impulse response of US system)
disp('Convolution:');
disp(conv(h,h));

a = 5; % Offset from edges
fprintf('a = %d\n',a)
fprintf('2+a = %d\n',2+a)
fprintf('2L-4-a = %d\n\n',2*L-4-a)

% Find (h*h)[2+a]
runTot = 0;
for k = 0:(2+a)
    if(2+a-k > L - 1 || 2+a-k < 0 || k > L -1 || k < 0)
        runTot = runTot +  0;
    else
        runTot = runTot + h(k + 1)*h(2+a-k+1);
    end
end
disp(sprintf('h*h[2+a]=%d',runTot));

% Find (h*h)[(2L-4) -a]
runTot = 0;
for k = 0:(2*L -4 -a) % %
    multPre = 1;
    arg1 = k;
    arg2 = 2*L-4-a-k;
    if(arg1> L - 1 || arg1 < 0 || arg2 > L -1 || arg2 < 0)
        runTot = runTot +  0;
    else
        runTot = runTot + multPre * h(arg1 + 1)*h(arg2 + 1);
    end
end
disp(sprintf('Old h*h[2L-4-a]=%d\n',runTot));

% Find (h*h)[(2L-4) -a] using new expression
runTot = 0;

if(L-1 < a-L+3)
    kIterSet = L-1:a-L+3;
else
    kIterSet = a-L+3:L-1;
end

for k = kIterSet % k = 0:(2*L -4 -a) % %
    multPre = 1;
    arg1 = k;
    arg2 = 2+a-k; %2*L-4-a-k; %2+a-k;
    if(arg1> L - 1 || arg1 < 0 || arg2 > L -1 || arg2 < 0)
        runTot = runTot +  0;
    else
        runTot = runTot + multPre * h(arg1 + 1)*h(arg2 + 1);
    end
    
end
disp(sprintf('New h*h[2L-4-a]=%d\n',runTot));

% Find (h*h)[(2L-4) -a] using newest expression
A = a + 2;
runTot = 0;
if(L-1 < A-(L-1))
    kIterSet = L-1:A-(L-1);
else
    kIterSet = A-(L-1):L-1;
end

for k = kIterSet % k = 0:(2*L -4 -a) % %
    multPre = 1;
    arg1 = k;
    arg2 = A-k; %2*L-4-a-k; %2+a-k;
    if(arg1> L - 1 || arg1 < 0 || arg2 > L -1 || arg2 < 0)
        runTot = runTot +  0;
    else
        runTot = runTot + multPre * h(arg1 + 1)*h(arg2 + 1);
    end
    
end
disp(sprintf('Newest h*h[2L-4-a]=%d\n',runTot));





