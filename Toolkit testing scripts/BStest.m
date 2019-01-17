% testing to check if using expmv works and to compare speeds

trunc = 20;
transmission = 0.5;    
theta = acos(sqrt(transmission));
z_abs = 0.1;
z_theta = 5.3;
alpha_abs = 0.13;
alpha_theta = 2.66;

psi_4 = TensorProduct({SingleModeSqueezedState([z_abs, z_theta], trunc), CoherentState([alpha_abs,alpha_theta],trunc), FockState(0,trunc)});
psi_4 = kron(eye(trunc),BeamSplitter(transmission,trunc)) * psi_4;

adag = TensorProduct({CreationOperator(trunc),eye(trunc),eye(trunc)});
a = TensorProduct({AnnihilationOperator(trunc),eye(trunc),eye(trunc)});
bdag = TensorProduct({eye(trunc^2), CreationOperator(trunc)});
b = TensorProduct({eye(trunc^2), AnnihilationOperator(trunc)});

A = - adag * b + bdag * a;

tic
BSMat = expm(theta * A);
BSMat = sparse(BSMat);

psi_4_original = BSMat * psi_4;
OriginalMethod = toc;

tic
psi_4_new = expmv(theta, A, psi_4,[],'double');
NewMethod = toc;

equaltol(psi_4_original, psi_4_new, 1e-4);