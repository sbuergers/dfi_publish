% David Meijer's customized function
% Based on PAL_PFML_negLLMultiple from Palamedes 1.8.2
%
% Computes (negative) Log Likelihood associated with 
% simultaneous fit of multiple Psychometric Functions
% using betabinomial model if desired

function negLL = DM_PFML_negLLMultiple_Beta(ParamsVector, paramsIDmatrix, StimLevels, NumPos, OutOfNum, PF, lapseLimits, betaBool)

%Transform the theta's (as Palamedes calls them) into a standard Palamedes parameter matrix [nCond x 4]     
if betaBool
    eta = ParamsVector(end);
    ParamsVector = ParamsVector(1:(end-1));
else
    eta = 0.5;  %This is a dummy, just to pass the parameter bounds check..
end
params = ParamsVector(paramsIDmatrix);

%Ensure parameter bounds 
if (min(params(:,4)) < lapseLimits(1)) || (max(params(:,4)) > lapseLimits(2)) || ...
        (eta >= 1) || (eta < 0) || ...
        (min(params(:,1)) < 0.025) || (max(params(:,1)) > 0.225) || ... % ensure threshold in measured range!
        (min(params(:,3)) < 0) || (max(params(:,3)) > 0.8) % ensure guess-rate between 0 and 0.8
    negLL = Inf;    
else
    %Compute Log-likelihood per condition
    nConditions = size(paramsIDmatrix,1);
    LLcond = zeros(1,nConditions);
    for c = 1:nConditions
        pcorrect = PF(params(c,:), StimLevels(c,:));
        if betaBool && (eta^2 > 10^-9) 
            %betabinomial model
            eta_prime = (1/eta^2)-1;
            alpha = pcorrect.*eta_prime;     
            beta = (1-pcorrect).*eta_prime;
            LLcond(c) = nansum(gammaln(NumPos(c,:)+alpha) ...
                       +gammaln((OutOfNum(c,:)-NumPos(c,:))+beta) ...
                       -gammaln(OutOfNum(c,:)+eta_prime) ...
                       +gammaln(eta_prime)-gammaln(alpha)-gammaln(beta));
        else
            %binomial model
            LLcond(c) = nansum(NumPos(c,:).*log(pcorrect)+(OutOfNum(c,:)-NumPos(c,:)).*log(1 - pcorrect));                  %Use nansum because according to MATLAB 0*log(0) = NaN.
        end
    end
    
    %Sum over conditions and change sign
    negLL = -sum(LLcond);
end

end %EOF