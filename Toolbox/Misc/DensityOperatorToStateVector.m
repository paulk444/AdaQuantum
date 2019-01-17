function state = DensityOperatorToStateVector (density_operator)

if not(PurityTest(density_operator, 0.01))
    disp('WARNING! Trying to output state vector when state is mixed!')
end

[V,D] = eig(full(density_operator));

[~,Loc] = ismembertol(1, real(diag(D)), 0.01);

state = V(:,Loc);

end