
trunc = randi(18) + 10;

alpha_abs = rand * 3;
alpha_phi = rand * 2 * pi;
z_abs = rand * 1.4;
z_phi = rand * 2 * pi;
phi = rand * 2 * pi;

psi_n = randi(10) - 1;
psi_alpha_abs = rand * 3;
psi_alpha_phi = rand * 2 * pi;
psi_z_abs = rand * 1.4;
psi_z_phi = rand * 2 * pi;

psi_in = TensorProduct( { FockState(psi_n, trunc), SingleModeSqueezedState([psi_z_abs, psi_z_phi], trunc), CoherentState([psi_alpha_abs, psi_alpha_phi], trunc) } );

f_old = @() OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, trunc);
f_new = @() NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, trunc);

TimeOld = timeit(f_old);
TimeNew = timeit(f_new);

results = [results; trunc, TimeOld, TimeNew];

%% Old method

function psi_out_old = OldMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, trunc)
A = DisplacementOperator([alpha_abs, alpha_phi], trunc);
B = PhaseShiftOperator(phi, trunc);
C = SingleModeSqueezingOperator([z_abs, z_phi], trunc);

psi_out_old = TensorProduct( { A, B, C } ) * psi_in;
end

%% New method

function psi_out_new = NewMethod(psi_in, alpha_abs, alpha_phi, z_abs, z_phi, phi, trunc)
psi_out_a = DisplacementOperation(psi_in, [alpha_abs, alpha_phi], trunc, 1, 3);
psi_out_b = PhaseShiftOperation(psi_out_a, phi, trunc, 2, 3);
psi_out_new = SingleModeSqueezingOperation(psi_out_b, [z_abs, z_phi], trunc, 3, 3);
end
