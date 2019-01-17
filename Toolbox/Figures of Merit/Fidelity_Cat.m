
function [Fidelity, parameters] = Fidelity_Cat(density_operator)

%======== Lewis O'Driscoll and Paul Knott, 2018 ============%
%========== https://arxiv.org/abs/1812.03183 ===============%

trunc = length(density_operator);

if trunc > 80
   disp("WARNING: Cat states in cat_states.mat were calculated with a truncation of 80"); 
end

% xs = linspace(-2, 2, 25);
% [xs, ys] = meshgrid(xs, xs);
% xs = reshape(xs, 1, numel(xs));
% ys = reshape(ys, 1, numel(ys));
% alphas = xs + 1i.*ys;
% 
% thetas = linspace(0, 2*pi, 10);
% [alphas, thetas] = meshgrid(alphas, thetas);
% 
% dims = size(alphas);
% fidlities = zeros(dims);
% for i = 1:dims(1)
%     for j = 1:dims(2)
%         alpha = alphas(i, j);
%         theta = thetas(i, j);
%         cat = CatState(alpha, theta, trunc);
%         fidelities(i, j) = abs(cat'*density_operator*cat);
%     end
% end

% Load search space of cat states from file
file = load('cat_states.mat');
states = file.cat_states;
alphas = file.alphas;
thetas = file.thetas;

% Calculate the fidelity with each cat state
fidelities = zeros(1, length(states));
for i = 1:length(states)
   cat = states(1:trunc, i);
   fidelities(i) = abs(cat'*density_operator*cat);
end

% Return the maximum fidelity
Fidelity = max(fidelities);

ind = find(fidelities == Fidelity, 1);
alpha = alphas(ind);
theta = thetas(ind);

parameters = cell({{'alpha', 'theta'}, {alpha, theta}});

end