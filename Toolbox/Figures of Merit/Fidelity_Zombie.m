
function [Fidelity, params] = Fidelity_Zombie(density_operator)

%======== Lewis O'Driscoll and Paul Knott, 2018 ============%
%========== https://arxiv.org/abs/1812.03183 ===============%

trunc = length(density_operator);

if trunc > 80
   disp("WARNING: States in zombie_states.mat were calculated with a truncation of 80"); 
end

% Load search space of cat states from file
file = load('zombie_states.mat');
states = file.zombie_states;
alphas = file.alphas;

% Calculate the fidelity with each cat state
fidelities = zeros(1, length(states));
for i = 1:length(states)
   zombie = full(states(1:trunc, i));
   fidelities(i) = abs(zombie'*density_operator*zombie);
end

% Return the maximum fidelity
Fidelity = max(fidelities);

ind = find(fidelities == Fidelity, 1);
alpha = alphas(ind);

params = cell({{'alpha'}, {alpha}});

end