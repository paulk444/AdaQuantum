function quadrature_eig = QuaderatureEigenstate(x_lambda, trunc)
    % Makes a quadrature eigenstate
    % Takes [|x_lambda|, lambda] as argument
    % Using Methods in Theoretical Quantum Optics by Barnett & Radmore A.4.11.
    
    r = x_lambda (1);
    lambda = x_lambda (2);
    
    a_dagger = CreationOperator(trunc);
    quadrature_eig = pi^(-1/4) * expm(-1/2 * r^2 * eye(trunc) + sqrt(2) * exp(1i*lambda) * r * a_dagger - 1/2 * exp(2i*lambda) * a_dagger^2) * FockState(0,trunc);
    
    
end
