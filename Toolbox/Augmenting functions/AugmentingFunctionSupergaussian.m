function AugmentingFunction = AugmentingFunctionSupergaussian (settings)
    
    steepness = settings(1);
    nbar_min = settings(2);
    nbar_max = settings(3);
    cutoff = settings(4);

    AugmentingFunction = @(nbar) exp(-2 * ((2 * nbar - nbar_min - nbar_max) / (nbar_max - nbar_min))^(2 * steepness)) * heaviside(nbar - cutoff);
    
end