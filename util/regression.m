function [r, b1, b0] = regression(x, y, varargin)
% Takes column vectors x and y as input and returns
% Pearson's correlation coefficient, slope and intercept
% parameters of a simple OLS regression.
%
% NOTE: This function is used as a proxy for the 'regression'
% function of the Neural Network Toolbox and therefore takes
% additional arguments that are unnecessary here.
    
    % remove NaNs
    keep = isfinite(y) & isfinite(x);
    x = x(keep);
    y = y(keep);
    
    % compute betas and pearson correlation
    p = polyfit(x, y, 1);
    b0 = p(2);  
    b1 = p(1);
    r = corr(x, y, 'type', 'Pearson', 'rows', 'complete');
    
end
