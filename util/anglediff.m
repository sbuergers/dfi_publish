function absDiffDeg = anglediff(a, b, dim)
    
    if dim == 1
        
        %Normalize the difference, abs operation is not necessary because mod(x,y) takes the sign of y.
        normDeg = mod(a-b,360);

        %This will be a number between 0-360, but we want the smallest angle which is between 0-180. Easiest way to get this is
        absDiffDeg = min(360-normDeg, normDeg);
        
    else
        
        error('Dim has to be 1...this function is just for johannas code, really...')
    
    end
    
return
