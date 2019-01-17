function projection_operator = HomodyneProjector(x_lambda, trunc)

    quadrature_eig = QuaderatureEigenstate(x_lambda, trunc);
    projection_operator = quadrature_eig * quadrature_eig';

end