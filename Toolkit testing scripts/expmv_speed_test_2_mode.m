
trunc = 120;

alpha_abs = rand * 3;
alpha_phi = rand * 2 * pi;
z_abs = rand * 1.4;
z_phi = rand * 2 * pi;

psi_alpha_abs = rand * 3;
psi_alpha_phi = rand * 2 * pi;
psi_z_abs = rand * 1.4;
psi_z_phi = rand * 2 * pi;

psi_in = TensorProduct( { SingleModeSqueezedState([psi_z_abs, psi_z_phi], trunc), CoherentState([psi_alpha_abs, psi_alpha_phi], trunc) } );

f_old = @() OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, trunc);
f_new = @() NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, trunc);

TimeOld = timeit(f_old);
TimeNew = timeit(f_new);

results = [results; trunc, TimeOld, TimeNew];

equaltol(feval(f_old),feval(f_new),1e-4);

%% Old method

function psi_out_old = OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, trunc)
A = DisplacementOperator([alpha_abs, alpha_phi], trunc);
B = SingleModeSqueezingOperator([z_abs, z_phi], trunc);

psi_out_old = TensorProduct( { A, B } ) * psi_in;
end

%% New method

function psi_out_new = NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, trunc)
psi_out_a = DisplacementOperation(psi_in, [alpha_abs, alpha_phi], trunc, 1, 2);
psi_out_new = SingleModeSqueezingOperation(psi_out_a, [z_abs, z_phi], trunc, 2, 2);
end
