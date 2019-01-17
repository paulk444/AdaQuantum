function new_results = expmv_speed_test_4_mode(results,trunc)


alpha_abs = rand * 3;
alpha_phi = rand * 2 * pi;
z_abs = rand * 1.4;
z_phi = rand * 2 * pi;
phi = rand * 2 * pi;
z_abs_2 = rand * 1.4;
z_phi_2 = rand * 2 * pi;

psi_n = randi(5) - 1;
psi_alpha_abs = rand * 3;
psi_alpha_phi = rand * 2 * pi;
psi_z_abs = rand * 1.4;
psi_z_phi = rand * 2 * pi;
psi_n_2 = randi(5) - 1;

psi_in = TensorProduct( { FockState(psi_n, trunc), SingleModeSqueezedState([psi_z_abs, psi_z_phi], trunc), CoherentState([psi_alpha_abs, psi_alpha_phi], trunc), FockState(psi_n_2, trunc) } );

f_old = @() OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, z_abs_2, z_phi_2, trunc);
f_new = @() NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, z_abs_2, z_phi_2, trunc);

% TimeOld = timeit(f_old);
TimeOld = 0;
TimeNew = timeit(f_new);

new_results = [results; trunc, TimeOld, TimeNew];

%equaltol(feval(f_old),feval(f_new),1e-4);

end

%% Old method

function psi_out_old = OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, z_abs_2, z_phi_2, trunc)
A = DisplacementOperator([alpha_abs, alpha_phi], trunc);
B = PhaseShiftOperator(phi, trunc);
C = SingleModeSqueezingOperator([z_abs, z_phi], trunc);
D = SingleModeSqueezingOperator([z_abs_2, z_phi_2], trunc);

psi_out_old = TensorProduct( { A, B, C, D } ) * psi_in;
end

%% New method

function psi_out_new = NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, z_abs_2, z_phi_2, trunc)
psi_out_a = DisplacementOperation(psi_in, [alpha_abs, alpha_phi], trunc, 1, 4);
psi_out_b = PhaseShiftOperation(psi_out_a, phi, trunc, 2, 4);
psi_out_c = SingleModeSqueezingOperation(psi_out_b, [z_abs, z_phi], trunc, 3, 4);
psi_out_new = SingleModeSqueezingOperation(psi_out_c, [z_abs_2, z_phi_2], trunc, 4, 4);
end
