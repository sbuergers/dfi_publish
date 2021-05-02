function SE_w = cousineau_within_subject_se(data)
% Compute Cousineau within subjec standard error to make
% figure error bars of within subject analyses more interpretable by eye.
% 
% INPUT
% data (matrix): The first dimension contains observations to compute SE
%     over, the second dimension denotes conditions and if present, the
%     third dimension denotes multiple cases for which to compute the ws
%     SE.

    subj_mean = nanmean(data,2);
    grand_mean = nanmean(subj_mean,1);
    data_corr = data - repmat(subj_mean, [1,size(data,2)]) + ...
        repmat(grand_mean, [size(data,1),size(data,2)]);
    
    SE_w = nanstd(data_corr,0,1)./sqrt(sum(~isnan(data_corr),1));
end



