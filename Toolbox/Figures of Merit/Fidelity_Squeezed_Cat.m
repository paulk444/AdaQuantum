function [Fidelity, params] = Fidelity_Squeezed_Cat(density_operator)

%======== Lewis O'Driscoll and Paul Knott, 2018 ============%
%========== https://arxiv.org/abs/1812.03183 ===============%

trunc = length(density_operator);

if trunc > 80
   disp("WARNING: States in squeezed_cat_states.mat were calculated with a truncation of 80"); 
end

% Load search space of cat states from file
file = load('squeezed_cat_states.mat');
states = file.sc_states;
alphas = file.alphas;
zs = file.zs;

% Calculate the fidelity with each cat state
fidelities = zeros(1, length(states));
for i = 1:length(states)
   sq_cat = full(states(1:trunc, i));
   fidelities(i) = abs(sq_cat'*density_operator*sq_cat);
end

% Return the maximum fidelity
Fidelity = max(fidelities);

ind = find(fidelities == Fidelity, 1);
alpha = alphas(ind);
z = zs(ind);

params = cell({{'alpha', 'z'}, {alpha, z}});

end