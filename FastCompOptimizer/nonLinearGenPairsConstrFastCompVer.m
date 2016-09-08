%clear
close all

% x holds all the codes
% Each two rows holds a pair

%intervals = [1 2 3 5 7 10 15];
intervals = [];

% f = The function we aim to minimize (maximize its reciprocal)
% Set to zero if you just want good autocorrelation
% Is ratio of cross correlation side lobe
% to autocorrelation main lobe
maxCC = @(x)maxXcorr(x, intervals);
mainLobe = @(x)minMainLobe(x, intervals);
fRatio =  @(x)maxCC(x)/mainLobe(x); 
fZero = @(x)0; % Use this if we don't care about cross correlation


%% Initial guess 

% Set code length
N = 5;

% Set number of pairs
numPairs = 3;

%startLoad = 1;
% Load in an appropriate number of pairs that satisfy
% the autocorrelation constraints
% createdLens = [10, 20, 50, 100];
% if (any(createdLens == N))    
%     fileName = strcat('C:\Users\Zemp-Lab\Desktop\OvernightPairGeneration\CompPairs_From_Nonlinear_Optimizer\\compPairs_len_',num2str(N),'.mat');
%     allPairsLoaded = load(fileName);
%     x0 = allPairsLoaded.('pairsSoFar');
% 
%     %startLoad = 23;
%     x0 = x0(startLoad:startLoad+2*numPairs -1,:);
% else
%     error('Codes of that length haven''t been made yet.');
% end
x0 = randn(numPairs*2, N);

% mu = 0; sigma = 1;  
% x0 = normrnd(mu,sigma,[2*numPairs,N]);  

%% Run algorithm

% No linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Bound constraints (how large can the values be in the codes?)
lb = -3*ones(size(x0(:)));
ub = 3*ones(size(x0(:)));

% Nonlinear constraint (ACF not at center goes to zero)
% Returns vector of ratios of sidelobes to mainlobes in ACfs
ACFConstr = @(x)ACFSumFuncConst(x,N,numPairs,intervals);

% A good display option is 'final'. Also 'iter'.
options = optimoptions('fmincon','Display','iter', 'UseParallel', 'always', 'MaxFunEvals',20000, 'algorithm', 'sqp');
%options = psoptimset('Display','iter', 'UseParallel', 'always');

% Set tolerance on constraint (not to be confused with CCF requirements)
% How big can the ratio of sidelobe / mainlobe in ACF be?
options.TolCon = 1e-4;

[x,~,~,output] = fmincon(fRatio,x0,A,b,Aeq,beq,lb, ub, ACFConstr, options);
metric = 1/fRatio(x);
%x = reshape(x, [5*2, 8]);

% If using NRI, add intervals of zeros between non-zero bits
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
    N = size(x, 2);
end

disp(['minAutoCorr/maxCrossCor = ', num2str(metric)]);
% Save the results
% disp(startLoad)
% save(strcat('NRI_Aug29',num2str(startLoad),'.mat'),'x')
save(strcat('len5_3codes.mat'),'x')

% What's the signal strength (main lobe) to interference (side lobe) ratio?
% if(numPairs > 1)
%     disp('Max main/side merit:')
%     disp(1/fRatio(x))
%     % Need to think this through
% %     disp('Welch best merit:')    
% %     disp(minMainLobe(x)/ welchBound( N, numPairs, 2 ))
% end

%% Plot results

plotRel = 0;

% Plot autocorrelation
if(plotRel == 1)
    % Create x axis (shift from -(codeLength-1) to +(codeLength-1)
    xAxis = -(N-1):(N-1);
    
    for i=1:numPairs
        currPair = x((2*i-1):2*i,:);
        plot(xAxis,xcorr(currPair(1,:)) + xcorr(currPair(2,:)));        
        hold on
    end    
    
    % Plot pairwise sum of cross correlation
    plotCcSum = 1;
    
    if (plotCcSum == 1 && numPairs > 1)
    
        pairChoices = nchoosek(1:numPairs,2);
        for i = 1:size(pairChoices,1)
            % We find the sum of pairwise cross correlations
            % Ex. for pair A and pair B, we sum xcorr(A1,B1) + xcorr(A2,B2)
            currChoices = pairChoices(i,:);
            firstPair = x((2*currChoices(1)-1):2*currChoices(1),:);
            secPair = x((2*currChoices(2)-1):2*currChoices(2),:);

            currXcorr = xcorr(firstPair(1,:),secPair(1,:)) + xcorr(firstPair(2,:),secPair(2,:));
            plot(xAxis,currXcorr)
            hold on
        end
        title(sprintf('ACF and CCF with Code Length = %d and #Pairs = %d', N,numPairs));
    else      
        title(sprintf('ACF with Code Length = %d and #Pairs = %d', N,numPairs));
    end
    
    xlabel('Shift amount')
    ylabel('Magnitude')    
end
