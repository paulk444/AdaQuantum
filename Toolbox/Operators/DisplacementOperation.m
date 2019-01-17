function  output_state = DisplacementOperation (input_state, alpha_abs_arg, trunc, acting_on_mode, number_of_modes)

% Acts on 'input_state' with a displacement operator, see https://arxiv.org/abs/1812.01032
% alpha_abs_arg(1) is the magnitude, alpha_abs_arg(2) is the phase
% This single-mode operator acts on mode 'acting_on_mode'
% The total number of modes must be specified. E.g. If we have a three mode system, and we want to act on mode 1,
% we choose: acting_on_mode = 1; number_of_modes = 3

%======== Rosanna Nichols and Paul Knott, 2018 =============%
%========== https://arxiv.org/abs/1812.01032 ===============%

    alpha = alpha_abs_arg(1) * exp(1i * alpha_abs_arg(2));
    exponent = alpha * CreationOperatorGeneral(trunc, acting_on_mode, number_of_modes) - alpha' * AnnihilationOperatorGeneral(trunc, acting_on_mode, number_of_modes);
    output_state = expmv(1, exponent, input_state, [], 'double');

end