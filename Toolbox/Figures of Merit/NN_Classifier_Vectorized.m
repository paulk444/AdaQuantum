function FOM_Values = NN_Classifier_Vectorized(state_vectors)
%NN_Classifier_Vectorised Converts an array of state vectors to an array of
%FOM values.
% Expects state_vectors to be a cell array, where state_vectors{n} is a
% 1x25 state vector

    [ver, exec, loaded] = pyversion;
    if (loaded == 0)
        pyversion '/usr/local/bin/python3'
    end

    net = py.quoptics.network.NeuralNetwork("model"); 
    
    % Convert MATLAB cell array -> numpy array
    state_vectors = py.list(state_vectors);
    state_vectors = py.numpy.array(state_vectors);
    
    preds = net.predict(state_vectors);
    preds = cell(preds);
    
    FOM_Values = zeros(1, length(state_vectors));
    for n = 1:length(preds)
       pred = preds{n};
       % Convert numpy array -> MATLAB array
       probs = pred{'probabilities'}.tolist();
       probs = cell2mat(cell(probs));
       % probs(6) is probability state is useless
       FOM_Values(n) = 1 - probs(6);
    end

end