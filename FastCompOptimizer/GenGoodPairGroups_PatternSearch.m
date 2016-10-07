clear
close all

% x holds all the codes
% Each two rows holds a pair

% f = The function we aim to minimize (maximize its reciprocal)
% Set to zero if you just want good autocorrelation
% Is ratio of cross correlation side lobe
% to autocorrelation main lobe
maxCC = @(x)maxXcorr(x);
mainLobe = @(x)minMainLobe(x);
fRatio =  @(x)maxCC(x)/mainLobe(x); 
fZero = @(x)0; % Use this if we don't care about cross correlation


%% Initial guess 
% Set code length
N = 100;

% Set number of pairs
numPairs = 2;

% Load in an appropriate number of pairs that satisfy
% the autocorrelation constraints
createdLens = [10, 20, 50, 100];
if (any(createdLens == N))    
    fileName = strcat('CompPairs_From_Nonlinear_Optimizer\\compPairs_len_',num2str(N),'.mat');
    allPairsLoaded = load(fileName);
    x0 = allPairsLoaded.('pairsSoFar');
    
    startLoad = 41;
    x0 = x0(startLoad:startLoad+2*numPairs,:);
else
    error('Codes of that length haven''t been made yet.');
end

% Randomly generate starting pairs
% mu = 0; sigma = 1;  
% x0 = normrnd(mu,sigma,[2*numPairs,N]);  

%%
% No linear constraints
A = []; b = []; Aeq = []; beq = []; 

% Bound constraints (how large can the values be in the codes?)
lb = []; % -3*ones(size(x0));
ub = []; %3*ones(size(x0));

% Nonlinear constraint (ACF not at center goes to zero)
% Returns vector of ratios of sidelobes to mainlobes in ACfs
ACFConstr = @(x)ACFSumFuncConst(x,N,numPairs);

maxNumIter = 2000000;
options = psoptimset('Display','iter',...
    'MaxFunEvals',maxNumIter, 'TolCon', 1e-4);

% Set tolerance on constraint (not to be confused with CCF requirements)
% How big can the ratio of sidelobe / mainlobe in ACF be?
[x,~,~,output] = patternsearch(fRatio,x0,A,b,Aeq,beq,lb, ub, ACFConstr, options);

% Save the results
% save('bestPairs.mat','x')

% What's the signal strength (main lobe) to interference (side lobe) ratio?
if(numPairs > 1)
    disp('Max main/side merit:')
    disp(1/fRatio(x))
    % Compare to Welch bound
    disp('Welch bound:')    
    disp(welchBound( N, numPairs, 2 ))
    
    % Would like to be welch bound
    % Max cross correlation of normalized codes
    % We set the magnitude of each individual pair to sqrt(2)
    disp('Welch metric:') 
    disp(maxCC(normr(x)/sqrt(2)))    
end

% Plot results
plotRel = 1;

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
