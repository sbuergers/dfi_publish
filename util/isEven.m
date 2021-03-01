function [ B ] = isEven( X )
%Takes integer as input and outputs a boolean depending on whether it is
%even or not

if mod(X,1) ~= 0
    B = false;
else

    if mod(X,2) == 0
        B = true;
    else
        B = false;
    end

end

end

