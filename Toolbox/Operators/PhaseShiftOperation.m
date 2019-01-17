function output_state = PhaseShiftOperation(input_state, phi, trunc, acting_on_mode, number_of_modes)

    output_state = expmv(phi, 1i * NumberOperator(trunc, acting_on_mode, number_of_modes), input_state, [], 'double');

end