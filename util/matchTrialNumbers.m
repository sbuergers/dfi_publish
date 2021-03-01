function [a, b] = matchTrialNumbers(a, b)
    if length(a) > length(b)
        a = randsample(a, length(b));
    else
        b = randsample(b, length(a));
    end    
end