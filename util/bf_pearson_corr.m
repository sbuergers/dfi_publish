function bf = bf_pearson_corr( r, n )
% Calculate bayes factors for correlation (compare M0: y=b0+e with M1:
% y=b0+b1*x+e
% Gives me exactly the same result as R, as well as the paper by Wetzel and
% Wagenmakers (2012), where the formula is from. So that's nice!
    fun = @(g) ( sqrt((n./2)) / gamma(1/2) )  .*  ...                          % term 1
               ( (1 + g).^((n-2)./2) )  .*  ...                                % term 2
               ( (1+(1-r.^2).*g).^-((n-1)./2).*g.^(-3/2).*exp((-n./(2.*g))) ); % term 3
    bf = integral(fun,0,Inf);
end












