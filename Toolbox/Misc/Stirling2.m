function S = Stirling2(n, k)
% Calculate Stirling number of the second kind

S = 0;
for j = 0:k
    S = S + (-1)^(k-j)*1/factorial(j)*1/factorial(k-j)*j^n;
end

end