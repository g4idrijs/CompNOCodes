clear
close all

bestSoFar = 0; % Best metric value so far

% Set NRI intervals
intervals = [];

% f = The function we aim to minimize (maximize its reciprocal)
% Set to zero if you just want good autocorrelation
% Is ratio of cross correlation side lobe
% to autocorrelation main lobe
maxCC = @(x)maxAllXCorr(x); % Consider cross correlation between all
mainLobe = @(x)minMainLobe(x,intervals);
fRatio =  @(x)maxCC(x)/mainLobe(x); 
fZero = @(x)0; % Use this if we don't care about cross correlation

%% Initial guess 
% Set code length
N = 100;

% Set number of pairs
numPairs = 2;

while(1 == 1)
    % Generate a set of numPairs complementary pairs with length N each
    x0 = [];
    for currPairInd = 1:numPairs
        x0 = [x0; genCompPair_opt(ones(1,N-1))];
    end
    x0 = x0/max(max(abs(x0)))*3;

    %%
    % No linear constraints
    A = []; b = []; Aeq = []; beq = []; 

    % Bound constraints (how large can the values be in the codes?)
    lb = []; % -3*ones(size(x0));
    ub = []; %3*ones(size(x0));

    % Nonlinear constraint (ACF not at center goes to zero)
    % Returns vector of ratios of sidelobes to mainlobes in ACfs
    ACFConstr = @(x)ACFSumFuncConst(x,N,numPairs,intervals);

    maxNumIter = inf;
    options = psoptimset('Display','iter',...
        'MaxFunEvals',maxNumIter, 'TolCon', 1e-4);

    % Set tolerance on constraint (not to be confused with CCF requirements)
    % How big can the ratio of sidelobe / mainlobe in ACF be?
    [x,~,~,output] = patternsearch(fRatio,x0,A,b,Aeq,beq,lb, ub, ACFConstr, options);
     metric = 1/fRatio(x);
    
    if(metric > bestSoFar)
        save(strcat('C:\Users\Zemp-Lab\Desktop\OvernightPairGeneration\GitCodes\CompNOCodes\Complementary Pairs\patternSearchAllXcorr\lowAllCC_2pairs_length100_',num2str(metric),'.mat'),'x')
        bestSoFar = metric;
    end
end
