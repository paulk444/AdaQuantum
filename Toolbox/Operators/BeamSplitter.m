function output_state = BeamSplitter(input_state, transmission, trunc, acting_on_modes, number_of_modes)

% Acts on 'input_state' with a beam splitter, see https://arxiv.org/abs/1812.01032
% This two-mode operator acts on modes 'acting_on_mode(1)' and 'acting_on_mode(2)'
% The total number of modes must be specified. E.g. If we have a three mode system, and we want to act with a beam splitter on modes 1
% and 2, we choose: acting_on_mode = [1 2]; number_of_modes = 3

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%


    theta = acos(sqrt(transmission));

    adag = CreationOperatorGeneral(trunc, acting_on_modes(1), number_of_modes);
    a = AnnihilationOperatorGeneral(trunc, acting_on_modes(1), number_of_modes);
    bdag = CreationOperatorGeneral(trunc, acting_on_modes(2), number_of_modes);
    b = AnnihilationOperatorGeneral(trunc, acting_on_modes(2), number_of_modes);
    
    exponent =  - adag * b + bdag * a;
    
    output_state = expmv(theta, exponent, input_state, [], 'double');
    
end