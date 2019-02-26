function B = flatten(A)
    %FLATTEN helper method to flatten matrix to 1D vector
    %flatten 2d-array to 1d
    AT = transpose(A); %transpose to flatten by row
    B = AT(:);
end

