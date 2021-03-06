function [SDM, SDM_adj] = dfi_SDM_2ifc( data, condh, condf , per_soa)
%
% takes dataset with 'acc' and 'trlid' columns as input, where trlid
% contains all relevant conditions for calculating perceptual sensitivity
% (dprime) and bias (C and lnB). 'condh' and 'condf' denote the values in
% 'trlid' to calculate hit and false alarm rate respectively. When
% 'per_soa' is false SOAs are disregarded in the calculations (default=1)
%
%            V2        V1
% resp '2'  P(H)      P(F)           P(H) = % V2 cor           condh
% resp '1'  P(M)      P(CR)          P(F) = % V1 inc           condf
%
%           V2A2      V1A2
% resp '2'  P(H)      P(F)           P(H) = % V2A2 cor         condh     9
% resp '1'  P(M)      P(CR)          P(F) = % V1A2 inc         condf     8
%
%           V1A1      V2A1
% resp '2'  P(M)      P(CR)          P(H) = % V1A1 cor         condh
% resp '1'  P(H)      P(F)           P(F) = % V2A1 inc         condf
%
% returns:
% 
%      pHit      pFA      d-prime   p Corr    crit C    crit lnB
%     0.6000    0.1000    1.5349    0.7500    0.5141    0.7891
%     0.8000    0.3000    1.3660    0.7500   -0.1586   -0.2167
%     0.9000    0.4000    1.5349    0.7500   -0.5141   -0.7891
%
% Stimulus overview
%     V A
% 1   0 0
% 2   1 0
% 3   2 0
% 4   0 1
% 5   1 1     % FUSION control
% 6   2 1     % FUSION  
% 7   0 2
% 8   1 2     % FISSION
% 9   2 2     % FISSION control
%10   1 2c (catch trials?)
%

if ~exist('per_soa', 'var');
    per_soa = 1;
end


if per_soa
    %% calculate per SOA: 
    %% pH, pF
    dall = data;
    soa  = unique(dall.soa);
    [ pF, pH, pF_adj, pH_adj ]  = deal(zeros(size(soa))); 
    
    % get hit and false alarm rates per soa
    for i = 1:numel(soa)
        d = dall(dall.soa == soa(i),:);
        pH(i) = sum(d.acc(ismember(d.trlid, condh)))    / sum(ismember(d.trlid, condh)); % perccor(V2A1)
        pF(i) = sum(d.acc(ismember(d.trlid, condf))==0) / sum(ismember(d.trlid, condf)); % percinc(V1A2)
        
        % Add 0.1 of each category (correct rejection, miss, hit, false
        % alarm)
        pH_adj(i) = (sum(d.acc(ismember(d.trlid, condh)))    + 0.1) / (sum(ismember(d.trlid, condh)) + 0.2);
        pF_adj(i) = (sum(d.acc(ismember(d.trlid, condf))==0) + 0.1) / (sum(ismember(d.trlid, condf)) + 0.2);
    end
    
    % dprime and criterion
    zH = PAL_PtoZ(pH);
    zF = PAL_PtoZ(pF);
    
    dP =  (zH - zF)./sqrt(2);
    C  = -(zH + zF)./sqrt(2);
    pC =  (pH + (1-pF))./2;
    
    % output 1
    SDM = dataset(soa, pH, pF, dP, pC, C, 'VarNames', {'soa', 'pH', 'pF', 'dP', 'pC', 'C'});
    
    % dprime and criterion, adjusted (avoid infinities)
    zH_adj = PAL_PtoZ(pH_adj);
    zF_adj = PAL_PtoZ(pF_adj);
    
    dP_adj =  (zH_adj - zF_adj)./sqrt(2);
    C_adj  = -(zH_adj + zF_adj)./sqrt(2);
    pC_adj =  (pH_adj + (1-pF_adj))./2;
    
    % output 2
    SDM_adj = dataset(soa, pH_adj, pF_adj, dP_adj, pC_adj, C_adj, 'VarNames', ...
        {'soa', 'pH_adj', 'pF_adj', 'dP_adj', 'pC_adj', 'C_adj'});
else
    error('Not yet implemented!')
end

end % funend

% eof









