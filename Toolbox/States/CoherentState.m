function [state] = CoherentState(alpha_abs_arg, trunc)
    % Creates coherent state. Takes [|alpha|, arg(alpha)] as argument

    alpha = alpha_abs_arg(1) * exp(1i * alpha_abs_arg(2));
    state = sparse(trunc,1);
    
    for n = 0:trunc-1
        state = state + alpha^n/sqrt(factorial(n)) * FockState(n, trunc);
    end                    
    
    state = state * exp(-abs(alpha)^2/2);
    
end