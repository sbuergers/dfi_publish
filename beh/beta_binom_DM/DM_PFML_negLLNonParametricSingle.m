% David Meijer's customized function
%
% Simple loop function:
% Calling PAL_PFML_negLLNonParametric multiple times,
% once for each condition (trial-type)

function [LL, nParams] = DM_PFML_negLLNonParametricSingle(NumPos, OutOfNum)

nTrialTypes = size(NumPos,1);

LL = 0;
nParams = 0;

for i=1:nTrialTypes
    [negLL, numParams] = PAL_PFML_negLLNonParametric(NumPos(i,:), OutOfNum(i,:));
    
    LL      = -1*negLL + LL;
    nParams = numParams + nParams;
end