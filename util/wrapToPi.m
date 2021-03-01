% Copied from the web. 
% lambdaWrapped = wrapToPi(lambda) wraps angles in lambda, in radians, to 
% the interval [?pi pi]. pi maps to pi and ?pi maps to ?pi. (In general, 
% odd, positive multiples of pi map to pi and odd, negative multiples of pi 
% map to ?pi.)
function xwrap = wrapToPi(x)

xwrap = rem (x, 2*pi);
idx = find (abs (xwrap) > pi);
xwrap(idx) = 2*pi * sign (xwrap(idx));


return