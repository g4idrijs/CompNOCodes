function genManyCodes()

% Generates complementary code pairs
% While minimizing the metric fRatio

% addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\Field_II')
% addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Field_II_Scripts\bft_64bit')
% field_init(-1);
% bft_init;

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
if(allAtOnce == 1) % Fire both in pair at same time
%     % Create image and measure metric
    fRatio = @(x)1/imageAndGetMetric(x); % Thing we want to minimize
    
% Minimizes cross correlation relative to autocorrelation between all pairs
    % Minimize cross correlation between codes in a pair
    % as well as the sum of cross correlation between pairs
%     if(isempty(intervals))
%         allAtOnceCC = @(x)maxAllXCorr(x);
%         fRatio =  @(x)(allAtOnceCC(x))/mainLobe(x);  
%     else
%         error('Intervals not supported for minimizing self cross correlation.')
%     end    
    
else % Multiple transmit events
    fRatio =  @(x)maxCC(x)/mainLobe(x); % Minimize cross correlation to autocorrelation between pairs
end

fZero = @(x)0; % Use this if we don't care about cross correlation

while(1 == 1)
    
    % Initial guess 
    % Set code length
    N = 3;
    % Set number of pairs that should get along
    numPairs = 1;
    
    % Random seed
    if(numPairs == 1)
        % Complementary pair
        addpath('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\genCompPairRecursive');
        x0 = genCompPair(ones(N-1,1)');   
    else
        % Random
        x0 = randn(numPairs*2, N);
    end    
    

    %% Run algorithm

    % No linear constraints
    A = []; b = []; Aeq = []; beq = []; 

    % Bound constraints (how large can the values be in the codes?)
    lb = -10*ones(size(x0(:)));
    ub = 10*ones(size(x0(:)));

    % Nonlinear constraint (ACF not at center goes to zero)
    % Returns vector of ratios of sidelobes to mainlobes in ACfs
    ACFConstr = @(x)ACFSumFuncConst(x,N,numPairs,intervals);

    % Get starting time (to keep track of how long we've been running)
    startTime = clock;    
    
    % A good display option is 'final'. Also 'iter'.
    options = optimoptions('fmincon','Display','iter', 'UseParallel', 'never',...
        'MaxFunEvals',5000000, 'algorithm', 'sqp', 'OutputFcn',@optTimer);
    %options = psoptimset('Display','iter', 'UseParallel', 'always');

    % Set tolerance on constraint (not to be confused with CCF requirements)
    % How big can the ratio of sidelobe / mainlobe in ACF be?
    options.TolCon = 5e-3; %1e-4
    
    % Optimize
    [x,~,~,output] = fmincon(fRatio,x0,A,b,Aeq,beq,lb, ub, ACFConstr, options);
%   
       

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
    
    % Save the results, labelled by the metric value
    metric = 1/fRatio(x); 
    if(metric > bestSoFar)
        save(strcat('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\FastCompOptimizer\allOnce_sigClutMetr\allOnce_sigClutMetr_',num2str(N)','_',num2str(metric),'.mat'),'x')
        bestSoFar = metric;
    end
end
    % Define a nested function to display information while the
    % optimization function is running (and save codes so far)
    function stop = optTimer(x,optimvalues,state)
        stop = false;
        
        % Display total elapsed time that it took to get to this iteration
        if isequal(state,'iter')
            disp(sprintf('Time so far: %.1f s',etime(clock,startTime)))
            % Store codes so far
            save(strcat('C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\FastCompOptimizer\allOnce_sigClutMetr\allOnce_sigClutMetr_','temp','.mat'),'x')

        end        

    end
end