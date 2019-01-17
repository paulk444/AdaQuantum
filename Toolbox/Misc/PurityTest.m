function isPure = PurityTest(density_operator, epsilon)
% Outputs bool: true if pure, false if not (with tolereance epsilon)
% Warning: this assumes the state is normalised!

    if abs(trace(density_operator^2)-1) < epsilon
       isPure = true;
    else
       isPure = false; 
    end

end