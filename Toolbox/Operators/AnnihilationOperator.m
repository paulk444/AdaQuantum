function [AnnihilationMatrix] = AnnihilationOperator(trunc)

% Creates the matrix form of the annihilation operator

i = zeros(trunc-1,1);
j = zeros(trunc-1,1);
v = zeros(trunc-1,1);
for k = 1:trunc-1
    i(k) = k;
    j(k) = k+1;
    v(k) = sqrt(k);
end

AnnihilationMatrix = sparse(i, j, v, trunc, trunc);

end