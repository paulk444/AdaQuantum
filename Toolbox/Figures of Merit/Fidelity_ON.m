function [Fidelity, params] = Fidelity_ON(density_operator)

%======== Lewis O'Driscoll and Paul Knott, 2018 ============%
%========== https://arxiv.org/abs/1812.03183 ===============%

% Find the truncation of the density operator
trunc = length(density_operator);
% Load the pre-calculated ON states
file = load('ON_states.mat');
states = file.ON_states;
ns = file.ns;
deltas = file.deltas;

fidelities = zeros(1, length(states));
% Calculate the fidelity of the density_operator with each of the ON states
for n = 1:length(states)
   state = states(1:trunc, n);
   fidelities(n) = abs(state'*density_operator*state);
end

% Return the maximum fidelity
Fidelity = max(fidelities);

ind = find(fidelities == Fidelity);
n = ns(ind);
delta = deltas(ind);

params = cell({{'N', 'delta'}, {n, delta}});

end
