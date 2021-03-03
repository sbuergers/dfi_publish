function dataPF = fitPSFmodel_sep_fits(data,betaBool)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute nr of right responses per location per StimType %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = data;

%SOAs used for this subject
SOAs = unique(d.soa)';
condition = unique(d.trlid);
NumCond   = length(condition);
numLevels = length(SOAs);
StimLevels = repmat(unique(d.soa)', [NumCond 1]);
OutOfNum = nan(NumCond,numLevels);            
NumCorr = nan(NumCond,numLevels);          
for icond=1:NumCond
    for j=1:numLevels
        accVect = d.acc(d.trlid == condition(icond) & d.soa == StimLevels(icond, j));
        NumCorr(icond, j) = sum(accVect);
        OutOfNum(icond, j) = length(accVect);
    end
end
nTrialsTotal = sum(sum(OutOfNum));
PropCorrect = NumCorr ./ OutOfNum;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare fitting parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Psychometric Function to fit
PF = @PAL_Weibull;

% Logistic (steeper than the cumulative normal)
sg.alpha   = 0.025:0.01:0.15;      % threshold (inflection point)
sg.beta    = 10.^[0.3010:0.25:3];  % slope (log scale)
sg.lambda  = 0.01:.015:0.15;        % lapse-rate
sg.gamma   = 0:.1:0.8;        % guess-rate
searchGrid = sg;
fparams    = [1 1 1 1];

% %Parameter grid defining parameter space through which to perform a brute-force
% %search for values to be used as initial guesses in iterative parameter search
% searchGrid.alpha  = linspace(-3,3,10);      %left-right bias (PSE)
% searchGrid.beta   = logspace(-1,1,10);      %slope (from 0.1 till 10, logarithmically spaced)
% searchGrid.gamma  = linspace(0,0.2,10);     %guess rate (this one needs to be defined, but since we use 'gammaEQlambda', it will not be used; searchGrid.lambda will be used instead) 
% searchGrid.lambda = linspace(0,0.2,10);     %lapse rate

%Hard limits for lapse rate (>0)
guessLimits = [0 0.8];                      %ensure guess-rate between 0 and 0.8
lapseLimits = [0 1];                        %lapse rate limits

%Nelder-Mead Simplex method options
searchOptions = PAL_minimize('options');    %type PAL_minimize('options','help') for help
searchOptions.MaxFunEvals = 5000;           %Allow loads of function evals (recommended minimum = 400*numParams(=9) = 3600)                      
searchOptions.MaxIter = 5000;               %Allow loads of iterations (there can be multiple (but at least one) function evaluations per iteration, so this maximum will never be reached)
searchOptions.TolFun = 1e-6;                %Increase desired precision on LL
searchOptions.TolX = 1e-6;                  %Desired precision on params (theta's)
searchOptions.Display = 'off';              %suppress fminsearch messages

%Number of bootstrap samples
BootstrapN = 1000;                          %Set to 5000 for accurate results             

%Create matrices: each row contains one condition
NumPos = NumCorr;

%For fits on individual conditions: Which parameters are free to vary (within the above constraints)?   
paramsFree = fparams;
% paramsFree = [1 1 1 1];                     %1: free parameter, 0: fixed parameter

%Fit independent PFs to each StimType in order to obtain initial parameter estimates for the dependent fits [below]    
paramsFitted_independent = nan(NumCond,4);
for i=1:NumCond
    paramsFitted_independent(i,:) = PAL_PFML_Fit(StimLevels(i,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,...
                                                'gammaEQlambda',0,'guessLimits',guessLimits,'lapseLimits',lapseLimits,'searchOptions',searchOptions);
end

%Specify the model and create a vector with initial parameter estimates
%10 parameter estimates per participant:
%3(cond, 2F, 2F1S, 1F2S) x 3(param, thresh., slope, guess-rate) + 1(param, lapse-rate)
%with ETA included, this actually makes 11 parameter estimates per participant!
paramsIDmatrix = 1:4;
InitialParamsVector = paramsFitted_independent;

%ENSURE INITIALIZATION IS WITHIN BOUNDS                 
lowerBounds = [0.025, -inf, guessLimits(1), lapseLimits(1)];
upperBounds = [0.025, inf, guessLimits(2), lapseLimits(2)];
InitialParamsVector = min(max(lowerBounds,InitialParamsVector),upperBounds);

%Betabinomial? 
if betaBool
    InitialParamsVector = [InitialParamsVector 0.05];                       %Initialize eta to 0.05
    lowerBounds = [lowerBounds 0];
    upperBounds = [upperBounds 1-1e-9];
end
nParams = length(InitialParamsVector);

%%%%%%%%%%%%%%%%%%%%%
%%% Fit the model %%%
%%%%%%%%%%%%%%%%%%%%%

%Perform model fit (dependent PF fits)
nFitTries = 0;
exitflag = 0;
InitialParamsVector2 = InitialParamsVector;
while (nFitTries < 10) && (exitflag == 0)
    [paramsFittedVector, negLL, exitflag, output] = ...
        PAL_minimize(@DM_PFML_negLLMultiple_Beta, InitialParamsVector2, searchOptions, ...
                     paramsIDmatrix, StimLevels, NumPos, OutOfNum, PF, lapseLimits, betaBool); 
    nFitTries = nFitTries+1;
    InitialParamsVector2 = min(max(lowerBounds,paramsFittedVector+1e-4),upperBounds);               %<--- ENSURE INITIALIZATION IS WITHIN BOUNDS
    
end %This while loop ensures that fminsearch converges. I found the tip here: https://uk.mathworks.com/matlabcentral/fileexchange/33328-improving-the-convergence-of-nelder-mead-and-so-fminsearch 

%Trow a warning if the fit did not converge
if ~(exitflag == 1)
    warning('MLE fit of model did not converge! Exitflag is %s', int2str(exitflag));
end

%Gather fitted parameters
if betaBool
    eta = paramsFittedVector(end);
    paramsFittedVector = paramsFittedVector(1:(end-1));
else
    eta = NaN;
end
paramsFitted = paramsFittedVector(paramsIDmatrix);
LL = -1*negLL;

%Compute likelihood of "saturated" model and deviance (for the Goodness-Of-Fit)
[LL_saturated, nParams_saturated] = DM_PFML_negLLNonParametricMultiple(NumPos, OutOfNum);
Dev = 2*(LL_saturated - LL);

%Compute theoretical Goodness-of-fit using the chi-square distribution
pDev_ChiSquare = 1-chi2cdf(Dev, nParams_saturated-nParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform a parametric bootstrap of the full model and check the Goodness-Of-Fit %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if BootstrapN > 0
    Bootstrap.NumPos       = nan(BootstrapN,NumCond,numLevels); 
    Bootstrap.paramsFitted = nan(BootstrapN,NumCond,4);                     
    Bootstrap.eta          = nan(BootstrapN,1);
    Bootstrap.LL           = nan(BootstrapN,1);                                 
    Bootstrap.Dev          = nan(BootstrapN,1); 
    Bootstrap.exitflag     = nan(BootstrapN,1); 
    for ib=1:BootstrapN

        %Simulate an experiment per StimType (using the parameters of the Full Model)  
        for i = 1:NumCond
            if betaBool
                Bootstrap.NumPos(ib,i,:) = DM_PF_SimulateObserverParametric_Beta(paramsFitted(i,:), StimLevels(i,:), OutOfNum(i,:), PF, eta);
            else
                Bootstrap.NumPos(ib,i,:) = PAL_PF_SimulateObserverParametric(paramsFitted(i,:), StimLevels(i,:), OutOfNum(i,:), PF);
            end
        end

        %Perform model fit (dependent PF fits)
        nFitTriesTmp = 0;
        exitflagTmp = 0;
        InitialParamsVector2 = InitialParamsVector;
        while (nFitTriesTmp < 10) && (exitflagTmp == 0)
            [paramsFittedVector, negLL, exitflagTmp, ~] = PAL_minimize(@DM_PFML_negLLMultiple_Beta, InitialParamsVector2, searchOptions, ...
                                                                      paramsIDmatrix, StimLevels, squeeze(Bootstrap.NumPos(ib,:,:))', ...
                                                                      OutOfNum, PF, lapseLimits, betaBool);  
            nFitTriesTmp = nFitTriesTmp+1;
            InitialParamsVector2 = min(max(lowerBounds,paramsFittedVector+1e-4),upperBounds);       %<--- ENSURE INITIALIZATION IS WITHIN BOUNDS
            
        end %This while loop ensures that fminsearch converges. I found the tip here: https://uk.mathworks.com/matlabcentral/fileexchange/33328-improving-the-convergence-of-nelder-mead-and-so-fminsearch

        %Gather fitted parameters
        if betaBool
            Bootstrap.eta(ib,1) = paramsFittedVector(end);
            paramsFittedVector = paramsFittedVector(1:(end-1));
        else
            Bootstrap.eta(ib,1) = NaN;
        end
        Bootstrap.paramsFitted(ib,:,:) = paramsFittedVector(paramsIDmatrix);
        Bootstrap.LL(ib,1) = -1*negLL;
        Bootstrap.exitflag(ib,1) = exitflagTmp;

        %Compute likelihood of "saturated" model and deviance (for the Goodness-Of-Fit)
        [LL_saturated, ~] = DM_PFML_negLLNonParametricMultiple(squeeze(Bootstrap.NumPos(ib,:,:))', OutOfNum);
        Bootstrap.Dev(ib,1) = 2*(LL_saturated - Bootstrap.LL(ib,1));
    end

    %Compute the Goodness-Of-Fit
    %Compare original Dev value with the bootstrapped distribution of Dev values (null hypothesis assumes that the PSF model is correct, since the bootstrap was based on those fitted parameters!)      
    %If we observe that only very few bootstrapped deviances are higher than the original (pDev < 0.05), than it is very unlikely that our original data was “created” using a model of the type that we fitted!    
    pDev_Bootstrap = (1+sum(Bootstrap.Dev > Dev)) / (1+BootstrapN);

    %Trow a warning if not all fits converged
    if sum(Bootstrap.exitflag) < BootstrapN
        warning('Only %s of %s Bootstrap model simulations converged!',int2str(sum(Bootstrap.exitflag)), int2str(BootstrapN));
    end
else
    Bootstrap = [];
    pDev_Bootstrap = [];
end % end of bootstrap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%%% End of function: Save data in struct %%%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 

%Collect all variables in struct 'dataPF' for the output of this function
dataPF.data.soa = SOAs;
dataPF.data.nTrialsTotal = nTrialsTotal;
dataPF.data.OutOfNum = OutOfNum;
dataPF.data.NumCorr = NumCorr;
dataPF.data.perCor = PropCorrect;

dataPF.settings.betaBool = betaBool;
dataPF.settings.BootstrapN = BootstrapN;
dataPF.settings.paramsIDmatrix = paramsIDmatrix;
dataPF.settings.nParams = nParams;
dataPF.settings.InitialParamsVector = InitialParamsVector;

dataPF.fit.nFitTries = nFitTries;
dataPF.fit.exitflag = exitflag;
dataPF.fit.output = output;
dataPF.fit.paramsFitted = paramsFitted;
dataPF.fit.eta = eta;
dataPF.fit.LL = LL;
dataPF.fit.Dev = Dev;
dataPF.fit.pDev_ChiSquare = pDev_ChiSquare;
dataPF.fit.pDev_Bootstrap = pDev_Bootstrap;

dataPF.Bootstrap = Bootstrap;

return %[EOF]
