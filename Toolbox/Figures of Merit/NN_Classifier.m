function output = NN_Classifier(density_operator)

[ver, exec, loaded] = pyversion;
if (loaded == 0)
    pyversion '/usr/local/bin/python3'
end

if not(PurityTest(density_operator, 0.01))
    % Can't do anything if we have a mixed state
    value = 0;
    return;
end

% Convert density matrix to a vector
state_vector = DensityOperatorToStateVector(density_operator);
% Network looks at absolute values of first 25 coefficients
state_data = zeros(1, 25);
if length(state_vector) >= 25
   state_data(1:25) = abs(state_vector(1:25));
   
else
    state_data(1:length(state_vector)) = abs(state_vector);
end
% Load the neural network
net = py.quoptics.network.NeuralNetwork("model");
% Classify the state
prob_dist = net.classify_dist(state_data);
% Convert python numpy array to MATLAB array
prob_dist = cell2mat(cell(prob_dist.tolist()));
% prob_dist(6) is probability state is useless - don't want this
output = max(prob_dist(1:5));

end