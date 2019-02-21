function B = flatten(A)
%flatten 2d-array to 1d
AT = transpose(A);
B = AT(:);
end

