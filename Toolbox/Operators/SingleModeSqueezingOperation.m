function output_state = SingleModeSqueezingOperation(input_state, z, trunc, acting_on_mode, number_of_modes)

% Acts on 'input_state' with a single-mode squeezing operator, see https://arxiv.org/abs/1812.01032
% z(1) is the magnitude of squeezing, z(2) is the squeezing phase
% This single-mode operator acts on mode 'acting_on_mode'
% The total number of modes must be specified. E.g. If we have a three mode system, and we want to act with a single-mode squeezing operator on mode 1
% , we choose: acting_on_mode = 1; number_of_modes = 3

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

    z = z(1) * exp(1i * z(2));

    exponent = - z/2 * CreationOperatorGeneral(trunc, acting_on_mode, number_of_modes)^2 + z'/2 * AnnihilationOperatorGeneral(trunc, acting_on_mode, number_of_modes)^2;
    output_state = expmv(1, exponent, input_state, [], 'double');

end