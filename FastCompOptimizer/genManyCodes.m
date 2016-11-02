% Generates complementary code pairs
% While minimizing the metric fRatio

clear 
bestSoFar = 0; % Best metric value so far

% f = The function we aim to minimize (maximize its reciprocal)
% Set to zero if you just want good autocorrelation
% Is ratio of cross correlation side lobe
% to autocorrelation main lobe
intervals = [];
maxCC = @(x)maxXcorr(x, intervals);
mainLobe = @(x)minMainLobe(x, intervals);

% Assume all codes are fired at once and minimize worse case cross
% correlation
allAtOnce = 1;
if(allAtOnce == 1)
    % Minimize cross correlation between codes in a pair
    % as well as the sum of cross correlation between pairs
    if(isempty(intervals))
        allAtOnceCC = @(x)maxAllXCorr(x);
        fRatio =  @(x)(allAtOnceCC(x))/mainLobe(x);  
    else
        error('Intervals not supported for minimizing self cross correlation.')
    end
else % Multiple transmit events
    fRatio =  @(x)maxCC(x)/mainLobe(x); % Minimize cross correlation to autocorrelation between pairs
end

fZero = @(x)0; % Use this if we don't care about cross correlation


while(1 == 1)
    
    % Initial guess 
    % Set code length
    N = 10;
    % Set number of pairs
    numPairs = 10;
    
    % Random seed
    x0 = randn(numPairs*2, N);

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
    options = optimoptions('fmincon','Display','iter', 'UseParallel', 'always', 'MaxFunEvals',5000000, 'algorithm', 'sqp');
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
    if(metric > bestSoFar)
        save(strcat('C:\Users\Zemp-Lab\Desktop\OvernightPairGeneration\GitCodes\CompNOCodes\FastCompOptimizer\lowAllCC_10pairs_length10\lowAllCC_10pairs_length10',num2str(metric),'.mat'),'x')
        bestSoFar = metric;
    end
end