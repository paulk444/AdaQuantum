function [state] = FockState(number, trunc)
                    
% Creates a Fock state with 'number' photons

    state = sparse(trunc,1);
    state(number+1) = 1;

                              
end