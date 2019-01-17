function [state] = TMFockState(number1,number2, trunc)

% Creates a two-mode Fock state with 'number1' photons in the first path, and 'number2' in the second

state = kron(FockState(number1,trunc),FockState(number2,trunc));

end