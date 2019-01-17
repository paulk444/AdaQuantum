function [CreationMatrix] = CreationOperator(trunc)

i = zeros(trunc-1,1);
j = zeros(trunc-1,1);
v = zeros(trunc-1,1);
for k = 1:trunc-1
    i(k) = k+1;
    j(k) = k;
    v(k) = sqrt(k);
end

CreationMatrix = sparse(i, j, v, trunc, trunc);

end