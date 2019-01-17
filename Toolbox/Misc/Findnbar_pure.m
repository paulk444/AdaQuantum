function nbar = Findnbar_pure(pure_state, number_of_modes)
% Takes state_vector as input, rather than a density matrix

    NumberOperator_full = sparse(length(pure_state),length(pure_state));
    
    trunc = round(length(pure_state)^(1/number_of_modes));
    
    for k = 1:number_of_modes
        
        NumberOperator_full = NumberOperator_full + NumberOperator(trunc, k, number_of_modes);
        
    end
    
    nbar = pure_state' * NumberOperator_full * pure_state;

end