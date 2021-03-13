% David Meijer's customized function
% Based on PAL_PF_SimulateObserverParametric from Palamedes 1.8.2
%
% Simulates responses using the parameters of the fitted PF. But every
% probability of the PF is drawn from a noisy beta distribution with the
% variance determined by "eta"

function NumPos = DM_PF_SimulateObserverParametric_Beta(paramsValues, StimLevels, OutOfNum, PF, eta)

%Gather probability of correct responses per stimulus level
pcorrect_raw = PF(paramsValues, StimLevels);

%Determine parameters of beta distribution for every stimulus level
eta_prime = (1/eta^2)-1;
alpha = pcorrect_raw.*eta_prime;                 
beta = (1-pcorrect_raw).*eta_prime;

%Draw the number of correct responses at random for each stimulus level
NumPos = zeros(1,length(StimLevels));
for iLevel = 1:length(StimLevels)
    
    %Draw pcorrect at random from beta distribution specific for this stimulus level
    pcorrect = betarnd(alpha(iLevel),beta(iLevel));
    
    %Draw random numbers from a uniform distribution between 0 and 1
    %The number of random numbers smaller than 'pcorrect' will be the simulated number of positive responses    
    Pos = rand(OutOfNum(iLevel),1);
    Pos(Pos < pcorrect) = 1;
    Pos(Pos ~= 1) = 0;
    NumPos(iLevel) = sum(Pos);
end