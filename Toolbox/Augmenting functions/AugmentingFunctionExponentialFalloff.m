function AugmentingFunction = AugmentingFunctionExponentialFalloff (settings)

    peak = settings(1);
    steepness = settings(2);
    cutoff = settings(3);

    M = exp(peak * steepness) * peak^(- peak * steepness);
    AugmentingFunction = @(nbar) M * exp(- steepness * nbar) * nbar^(peak * steepness) * heaviside(nbar - cutoff);
    
end