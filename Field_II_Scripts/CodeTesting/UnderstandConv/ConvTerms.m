clear
% Impulse response of transducer
L = 6;
n = 0:L-1;
m = 2;
h = sin(2*pi*m/(L-1)*n);
% h = [0 1 2 0]; 
L = numel(h);

% Self convolution (proportional to impulse response of US system)
% disp('Convolution:');
% disp(conv(h,h));

a = 2*6-4; % Offset from edges
A = a + 2;
fprintf('a = %d\n',a)
fprintf('2+a = %d\n',2+a)
fprintf('2L-4-a = %d\n\n',2*L-4-a)

% Find (h*h)[A]
currTermsOld = [];
runTot = 0;
for k = 0:A
    if(A-k > L - 1 || A-k < 0 || k > L -1 || k < 0)
        runTot = runTot +  0;
        currTerm = 0;
    else
        runTot = runTot + h(k + 1)*h(A-k+1);
        currTerm = h(k + 1)*h(A-k+1);
    end
    currTermsOld = [currTermsOld, currTerm];
end
disp(sprintf('h*h[2+a]=%d',runTot));
[currTermsOld', (0:A)']

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


% Find (h*h)[(2L-4) -a] using newest expression
runTot = 0;
if(L-1 < A-(L-1))
    kIterSet = L-1:A-(L-1);
else
    kIterSet = A-(L-1):L-1;
end

currTermsNew = [];
for k = kIterSet 
    multPre = 1;
    arg1 = k;
    arg2 = A-k; %
    if(arg1> L - 1 || arg1 < 0 || arg2 > L -1 || arg2 < 0)
        runTot = runTot +  0;
        currTerm = 0;
    else
        runTot = runTot + multPre * h(arg1 + 1)*h(arg2 + 1);
        currTerm = multPre * h(arg1 + 1)*h(arg2 + 1);
    end
    currTermsNew = [currTermsNew currTerm];
end
disp(sprintf('Newest h*h[2L-4-a]=%d\n',runTot));

[currTermsNew', (kIterSet)']

% Find (h*h)[A] using universal expression
runTot = 0;
kIterSet = max(1,A-L+2):min(A,L-2);

currTermsUniv= [];
for k = kIterSet 
    multPre = 1;
    arg1 = k;
    arg2 = A-k; %

    runTot = runTot + multPre * h(arg1 + 1)*h(arg2 + 1);
    currTerm = multPre * h(arg1 + 1)*h(arg2 + 1);
   
    currTermsUniv = [currTermsUniv currTerm];
end
disp(sprintf('Universal h*h[A]=%d\n',runTot));

[currTermsUniv', kIterSet']





