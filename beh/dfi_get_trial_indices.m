function [ id, idm, idconf, idconfm ] = dfi_get_trial_indices( beh, conf_judg )
%[ id, idm, idconf, idconfm ] = dfi_get_trial_indices( beh, conf_judg )
% If conf_judg is 'yes' instead of 'no' there are 8 response buttons
% instead of two.
%   load trial indices of dfi data set 'beh'
%       
%    STIM_RESP
%         v1: [0x1 double]
%      v1_F2: [0x1 double]
%      v1_F1: [0x1 double]
%         v2: [0x1 double]
%      v2_F2: [0x1 double]
%      v2_F1: [0x1 double]
%         a1: [0x1 double]
%      a1_S2: [0x1 double]
%      a1_S1: [0x1 double]
%         a2: [0x1 double]
%      a2_S2: [0x1 double]
%      a2_S1: [0x1 double]
%        av1: [84x1 double]
%     av1_F2: [8x1 double]
%     av1_F1: [76x1 double]
%        fus: [82x1 double]
%     fus_F2: [26x1 double]
%     fus_F1: [56x1 double]
%        fis: [84x1 double]
%     fis_F2: [39x1 double]
%     fis_F1: [45x1 double]
%        av2: [84x1 double]
%     av2_F2: [45x1 double]
%     av2_F1: [39x1 double]

%% 2.) --- TRIAL INDICES ---
if ~exist('conf_judg', 'var')
    conf_judg = 'no';
end

% Make sure responses are adjusted for response keys
if strcmp(conf_judg, 'no')
    resp = beh.resp;
    resp(beh.resp == 1 & strcmp(beh.rkeycfg, 'BA')) = 2;
    resp(beh.resp == 2 & strcmp(beh.rkeycfg, 'BA')) = 1;
    beh.resp = resp;
elseif strcmp(conf_judg, 'yes')
    resp = beh.resp;
    resp(beh.resp <= 4 & strcmp(beh.rkeycfg, 'AB')) = 1;
    resp(beh.resp >= 5 & strcmp(beh.rkeycfg, 'AB')) = 2;
    resp(beh.resp <= 4 & strcmp(beh.rkeycfg, 'BA')) = 2;
    resp(beh.resp >= 5 & strcmp(beh.rkeycfg, 'BA')) = 1;
    beh.resp = resp;
end

% unimodal trials
id.v1     = find(beh.trlid == 2);                      % 1 Flash
id.v1_F2  = find(beh.trlid == 2 & beh.resp == 2);      % 1 Flash report 2
id.v1_F1  = find(beh.trlid == 2 & beh.resp == 1);      % 1 Flash report 1
id.v2     = find(beh.trlid == 3);
id.v2_F2  = find(beh.trlid == 3 & beh.resp == 2);
id.v2_F1  = find(beh.trlid == 3 & beh.resp == 1);
id.a1     = find(beh.trlid == 4);
id.a1_S2  = find(beh.trlid == 4 & beh.resp == 2);
id.a1_S1  = find(beh.trlid == 4 & beh.resp == 1);
id.a2     = find(beh.trlid == 7);
id.a2_S2  = find(beh.trlid == 7 & beh.resp == 2);
id.a2_S1  = find(beh.trlid == 7 & beh.resp == 1);
% multimodal trials
id.av1    = find(beh.trlid == 5);
id.av1_F2 = find(beh.trlid == 5 & beh.resp == 2);
id.av1_F1 = find(beh.trlid == 5 & beh.resp == 1);
id.fus    = find(beh.trlid == 6);
id.fus_F2 = find(beh.trlid == 6 & beh.resp == 2);
id.fus_F1 = find(beh.trlid == 6 & beh.resp == 1);
id.fis    = find(beh.trlid == 8);
id.fis_F2 = find(beh.trlid == 8 & beh.resp == 2);
id.fis_F1 = find(beh.trlid == 8 & beh.resp == 1);
id.av2    = find(beh.trlid == 9);
id.av2_F2 = find(beh.trlid == 9 & beh.resp == 2);
id.av2_F1 = find(beh.trlid == 9 & beh.resp == 1);


% match trial numbers in second structure
idm = id;
[idm.v1_F2, idm.v1_F1]   = matchTrialNumbers(idm.v1_F2, idm.v1_F1);
[idm.av1_F2, idm.av1_F1] = matchTrialNumbers(idm.av1_F2, idm.av1_F1);
[idm.v2_F2, idm.v2_F1]   = matchTrialNumbers(idm.v2_F2, idm.v2_F1);
[idm.fis_F2, idm.fis_F1] = matchTrialNumbers(idm.fis_F2, idm.fis_F1);
[idm.fus_F2, idm.fus_F1] = matchTrialNumbers(idm.fus_F2, idm.fus_F1);
[idm.av2_F2, idm.av2_F1] = matchTrialNumbers(idm.av2_F2, idm.av2_F1);


% add confidence
idconf = [];
if strcmp(conf_judg, 'yes')
    % unimodal trials
    idconf.v1_12  = find(beh.trlid == 2 & ismember(beh.conf, [1,2]));  
    idconf.v1_34  = find(beh.trlid == 2 & ismember(beh.conf, [3,4]));
    idconf.v2_12  = find(beh.trlid == 3 & ismember(beh.conf, [1,2]));
    idconf.v2_34  = find(beh.trlid == 3 & ismember(beh.conf, [3,4]));
    idconf.a1_12  = find(beh.trlid == 4 & ismember(beh.conf, [1,2]));
    idconf.a1_34  = find(beh.trlid == 4 & ismember(beh.conf, [3,4]));
    idconf.a2_12  = find(beh.trlid == 7 & ismember(beh.conf, [1,2]));
    idconf.a2_34  = find(beh.trlid == 7 & ismember(beh.conf, [3,4]));
    % multimodal trials
    idconf.av1_12 = find(beh.trlid == 5 & ismember(beh.conf, [1,2]));
    idconf.av1_34 = find(beh.trlid == 5 & ismember(beh.conf, [3,4]));
    idconf.fus_12 = find(beh.trlid == 6 & ismember(beh.conf, [1,2]));
    idconf.fus_34 = find(beh.trlid == 6 & ismember(beh.conf, [3,4]));
    idconf.fis_12 = find(beh.trlid == 8 & ismember(beh.conf, [1,2]));
    idconf.fis_34 = find(beh.trlid == 8 & ismember(beh.conf, [3,4]));
    idconf.av2_12 = find(beh.trlid == 9 & ismember(beh.conf, [1,2]));
    idconf.av2_34 = find(beh.trlid == 9 & ismember(beh.conf, [3,4]));
    
    % match trial numbers
    idconfm = idconf;
    [idconfm.v2_12, idconfm.v2_34]   = matchTrialNumbers(idconfm.v2_12, idconfm.v2_34);
    [idconfm.fis_12, idconfm.fis_34] = matchTrialNumbers(idconfm.fis_12, idconfm.fis_34);
    [idconfm.fus_12, idconfm.fus_34] = matchTrialNumbers(idconfm.fus_12, idconfm.fus_34);
    [idconfm.av2_12, idconfm.av2_34] = matchTrialNumbers(idconfm.av2_12, idconfm.av2_34);
end


% NESTED FUNCTION
function [a, b] = matchTrialNumbers(a, b)
    if length(a) > length(b)
        a = randsample(a, length(b));
    else
        b = randsample(b, length(a));
    end    
end


end



% Giulio!

% see1 = find(((dall.resp==1 & strcmp(dall.rkeycfg,'AB')) | (dall.resp==2 & strcmp(dall.rkeycfg,'BA'))) & dall.trlid == 3);
% see2 = find(((dall.resp==2 & strcmp(dall.rkeycfg,'AB')) | (dall.resp==1 & strcmp(dall.rkeycfg,'BA'))) & dall.trlid == 3);



