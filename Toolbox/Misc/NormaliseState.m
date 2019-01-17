function normalised_state = NormaliseState(unnormalised_state)

% Takes pure state and normalises it

normalisation_constant = (unnormalised_state' * unnormalised_state) ^ (-1/2);

normalised_state = normalisation_constant * unnormalised_state;