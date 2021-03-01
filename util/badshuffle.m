function nogoodshuffle = badshuffle(mat, nrep)
% goes through a matrix of numerical items and checks if there are a number
% of nrep repititions in succession going from row to row. For this all
% elements within one row are summed. If there are it gives 1, otherwise 0.
csum = 1;
msum = 1;
vect = sum(mat,2);
if numel(vect) <= nrep
    nogoodshuffle = 0;
    return; 
end;
for i = 2:numel(vect)
    if vect(i) == vect(i-1)
        csum = csum + 1;
        if csum > msum
            msum = csum;
        end
    else
        csum = 1;
    end
end % efor
if msum >= nrep
    nogoodshuffle = 1;
else
    nogoodshuffle = 0;
end
end % efun