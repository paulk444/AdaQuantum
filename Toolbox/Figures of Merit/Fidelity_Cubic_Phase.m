
function [Fidelity, params] = Fidelity_Cubic_Phase(density_operator)

%======== Lewis O'Driscoll and Paul Knott, 2018 ============%
%========== https://arxiv.org/abs/1812.03183 ===============%

trunc = length(density_operator);

if trunc > 80
   disp("WARNING: States in cubic_phase_states.mat were calculated with a truncation of 80"); 
end

% Load search space of states from file
file = load('cp_states.mat');
states = file.cp_states;
gammas = file.gammas;
rs = file.rs;

% Restrict search to gamma >= 0.05
states = states(:, gammas >= 0.05);
rs = rs(gammas >= 0.05);
gammas = gammas(gammas >= 0.05);

% Calculate the fidelity with each state
fidelities = zeros(1, length(states));
for i = 1:length(states)
   cubic_phase = full(states(1:trunc, i));
   fidelities(i) = abs(cubic_phase'*density_operator*cubic_phase);
end

% Return the maximum fidelity
Fidelity = max(fidelities);

ind = find(fidelities == Fidelity, 1);
gamma = gammas(ind);
r = rs(ind);

params = cell({{'gamma', 'r'}, {gamma, r}});

end