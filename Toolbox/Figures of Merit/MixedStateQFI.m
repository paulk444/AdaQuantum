function [fisherinfo] = MixedStateQFI(density_operator)

% Calculates the mixed state QFI using eqn (4.18) in Paul Knott's thesis: http://etheses.whiterose.ac.uk/8931/1/Paul%20Thesis%20Final%202015.pdf

%=========== Lana Mineh and Paul Knott 2018 ================%
%========== https://arxiv.org/abs/1812.01032 ===============%

trunc = length(density_operator);
N = NumberOperator(trunc, 1, 1);

% Calculate eigenvalues and eigenvectors
[eigenvectors, eigenvalues] = eig(full(density_operator));
eigenvalues = real(diag(eigenvalues));


% Test eigenvectors are orthogonal. If not, make rho Hermitian

sum_orth=0;

for i=1:trunc
    for j=1:trunc
        sum_orth=sum_orth+eigenvectors(:,i)'*eigenvectors(:,j);
    end
end

if abs(sum_orth-trunc) > 10^(-6)
    %fprintf('Eigenvectors not initially orthogonal \n')
    density_operator=(density_operator+density_operator')/2;
    [eigenvectors, eigenvalues] = eig(full(density_operator));
    eigenvalues = real(diag(eigenvalues));
end


% Calculate QFI
fisherinfo = 0;
epsilon = 1e-6;
for i = 1:trunc
    if eigenvalues(i) > epsilon
        fisherinfo = fisherinfo + eigenvalues(i)*PureStateQFI(eigenvectors(:,i)*eigenvectors(:,i)');
        
        for j = 1:trunc
            if (i ~= j) && (eigenvalues(j) > epsilon)
                element = eigenvectors(:,i)'*N*eigenvectors(:,j);
                fisherinfo = fisherinfo - 8*eigenvalues(i)*eigenvalues(j)/(eigenvalues(i) + eigenvalues(j))*(element*conj(element));
            end
        end
    end
end

fisherinfo = real(fisherinfo);

end




