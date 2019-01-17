function output_state = ApplyLossBS(input_state, loss_rate)

% Function to apply loss using a beam splitter to a one mode pure state

%=========== Lana Mineh and Paul Knott 2018 ================%
%========== https://arxiv.org/abs/1812.01032 ===============%


trunc = length(input_state);
total_input = kron(input_state, FockState(0,trunc));
Id = speye(trunc);

% Apply beam-splitter
outBS = BeamSplitter(1-loss_rate, trunc)*total_input;
rho = outBS*outBS';

% Trace over mode 2 to get output density matrix
output_state = sparse(trunc,trunc);
for i = 0:trunc-1
    Fock_i = kron(Id, FockState(i,trunc));
    output_state = output_state + Fock_i'*rho*Fock_i;
end

end