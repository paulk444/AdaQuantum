function nbar = Findnbar(density_operator, number_of_modes)

    NumberOperator_full = sparse(length(density_operator),length(density_operator));
    
    trunc = round(length(density_operator)^(1/number_of_modes));
    
    for k = 1:number_of_modes
        
        NumberOperator_full = NumberOperator_full + NumberOperator(trunc, k, number_of_modes);
        
    end
    
    nbar = trace(density_operator * NumberOperator_full);

end